#include <spdlog/spdlog.h>
#include "iostream"

#include "scene.hpp"
#include "IO/MeshIO.hpp"
#include "Io/FileUtils.hpp"

#include "Bsdfs/Reflection.hpp"
#include "Bsdfs/BsdfFactory.hpp"

#include "Primitives/TriangleMesh.hpp"
#include "Primitives/Sphere.hpp"
#include "Primitives/Quad.hpp"
#include "Primitives/Cube.hpp"

#include "Lights/Infinte.hpp"

static int occludedCount = 0;
static int occludedCount1 = 0;


Scene::Scene(const nlohmann::json j) {
//    std::unordered_map<std::string,nlohmann::json> bsdf_jsons =j.at("bsdfs");
//
//    std::vector<nlohmann::json> bsdf_jsons = j.at("bsdfs");
//
//    std::unordered_map<std::string,std::shared_ptr<Bsdf>> bsdfs;
    bsdfs = BSDFFactory::LoadBsdfsFromJson(j.at("bsdfs"));
//    for(auto item:bsdf_jsons){
//        bsdfs[item.first]=BsdfFactory::LoadBsdfFromJson(item.second);
//    }

    auto loadSphere = [](const nlohmann::json & json, std::shared_ptr < BSDF > bsdf) {
        return std::make_shared < Sphere >(json.at("radius"), bsdf);
    };

    auto loadQuad = [](const nlohmann::json & json, std::shared_ptr < BSDF > bsdf) {
        return std::make_shared < Quad >(bsdf);
    };
    auto loadCube = [](const nlohmann::json & j, std::shared_ptr < BSDF > bsdf) {
        return std::make_shared < Cube >(bsdf);
    };
    auto loadInfiniteSphere = [](const nlohmann::json & j, std::shared_ptr < BSDF > bsdf) {
        return nullptr;
    };


    std::unordered_map < std::string,
            std::function < std::shared_ptr < Primitive >(const nlohmann::json, std::shared_ptr < BSDF > bsdf) >
    > loadMap{{"sphere",          loadSphere},
              {"quad",            loadQuad},
              {"cube",            loadCube},
              {"infinite_sphere", loadInfiniteSphere},
    };


    for ( auto p: j.at("primitives") ) {
        //get bsdf
        std::shared_ptr < BSDF > bsdf = fetchBSDFFromJson(p.at("bsdf"));


        Transform * transform = nullptr;
        if ( p.contains("transform") )
            transform = new Transform(p["transform"]);


        std::string type = p.at("type");
        spdlog::info(type);
        if ( type == "mesh" ) {
            std::shared_ptr < TriangleMesh > mesh = std::make_shared < TriangleMesh >();
            mesh->Load(p, * this, transform);
            primitives.push_back(mesh);
            // primitives.insert(primitives.end(),triangles.begin(),triangles.end());
            handleAddLight(p, primitives.size() - 1, primitives.size());
        } else {
            std::shared_ptr < Primitive > primitive = loadMap[type](p, bsdf);
            if ( transform && primitive ) primitive->transform(* transform);
            if ( primitive && type != "cube"
           // && (primitive->bsdf->name=="backWall" ||primitive->bsdf->name=="light" )
           // && type!="quad"
            ) {
                primitives.push_back(primitive);
                handleAddLight(p, primitives.size() - 1, primitives.size());
            } else {
                if ( type == "infinite_sphere" ) {
                    auto bitMap = std::make_shared < BitMapTexture >
                            (FileUtils::WorkingDir + p.at("emission").get < std::string >());
                    lights.push_back(std::make_shared < InfinteSphere >(bitMap));
                }
            }
        }
    }

//    primitives=std::vector<std::shared_ptr<Primitive>>(primitives.begin(),primitives.begin()+2);

    spdlog::info("{} Primitives", primitives.size());
    spdlog::info("{} lights", lights.size());
    spdlog::info("{} Bsdfs", bsdfs.size());

    _useBVH = getOptional(j["renderer"], "scene_bvh", true);
   // _useBVH = false;

}


void Scene::handleAddLight(const nlohmann::json & p, size_t l, size_t r) {
    if ( p.find("light") != p.end() ) {
        auto lightJson = p.at("light");
        std::string lightType = lightJson.at("type");

        if ( lightType == "area" ) {
            Spectrum albedo = lightJson.at("albedo");
            Float totalArea = 0, invArea;
            for ( size_t i = l ; i < r ; i ++ ) totalArea += primitives[i]->Area();
            invArea = 1 / totalArea;
            for ( size_t i = l ; i < r ; i ++ ) {
                // flux to radiosity
                auto light = std::make_shared < AreaLight >(this->primitives[i], albedo * invArea);

                this->primitives[i]->areaLight = light;
                lights.push_back(light);
            }
        }
    }
    Spectrum emission;
    if ( containsAndGet(p, "emission", emission) ) {
        for ( size_t i = l ; i < r ; i ++ ) {
            auto light = std::make_shared < AreaLight >(this->primitives[i], emission);
            this->primitives[i]->areaLight = light;
            lights.push_back(light);
        }
    }
    return;
}


std::optional < Intersection > Scene::intersect(const Ray & ray) const {
    if ( _useBVH ) {
        Ray tempRay(ray);
        RTCRayHit rayHit;
        EmbreeUtils::convertRay(& ray, & rayHit);
        //return bvh->intersect(ray);
        RTCIntersectContext context;
        rtcInitIntersectContext(& context);
        rtcIntersect1(_scene, &context, &rayHit);
        if(ray.farT == 100)
            int k=1;
        if ( rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID ) {
            return std::nullopt;
        }
        if(ray.farT == 100)
            int k=1;
        Intersection * its = EmbreeUtils::RTCRayHit1_(& rayHit)->its;
        its->w =  ray.d;
        return {*its};
    }

    std::optional < Intersection > minIntersection;
    Ray _ray(ray);
    for ( auto primitive: primitives ) {

        auto its = primitive->intersect(_ray);
        if ( its.has_value() ) {
            minIntersection = its;
        }
    }
    if ( minIntersection.has_value() ) {
        minIntersection->w = ray.d;
    }
    return minIntersection;
}

bool Scene::intersectP(const Ray & ray) const {
    if ( _useBVH ) {
        RTCRay shadowRay;
        EmbreeUtils::convertRay(& ray, & shadowRay);
        RTCIntersectContext context;
        rtcInitIntersectContext(& context);
        rtcOccluded1(_scene, &context, &shadowRay);
        return  shadowRay.tfar == _NEG_INFINY ;
    } else {
        Ray _ray(ray);
        for ( auto primitive: primitives ) {
            if ( primitive->occluded(ray) )
                return true;
        }
    }

    return false;
}

void Scene::logDebugInfo( ) {
    spdlog::info("Occluded {0} {1}", occludedCount, occludedCount1);
//    static_cast<InfinteSphere *>(lights[0].get())->logDebugInfo();
}


void Scene::setUp( ) {


    if ( _useBVH ) {
        _scene = rtcNewScene(EmbreeUtils::getDevice());
        auto temp = rtcGetDeviceError(EmbreeUtils::getDevice());
        for ( const auto & primitve: primitives ) {
            RTCGeometry geometry = primitve->initRTC();
            auto geom = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_USER);
            rtcEnableGeometry(geom);
            rtcSetGeometryUserPrimitiveCount(geom, 1);
            rtcSetGeometryUserData(geom, primitve.get());
            rtcSetGeometryBoundsFunction(geom, & EmbreeUtils::instanceBoundsFunc, nullptr);
            auto error = rtcGetDeviceError(EmbreeUtils::getDevice());
            rtcSetGeometryIntersectFunction(geom, & EmbreeUtils::instanceIntersectFunc);
            error = rtcGetDeviceError(EmbreeUtils::getDevice());
            rtcSetGeometryOccludedFunction(geom, & EmbreeUtils::instanceOccludedFunc);
            error = rtcGetDeviceError(EmbreeUtils::getDevice());
            rtcCommitGeometry(geom);
            int res = rtcAttachGeometry(_scene, geom);
            spdlog::info(res);
        }
        rtcCommitScene(_scene);
        //bvh = std::make_unique<BVHAccel>(primitives);  // I will use embree now ,the bvh I implmented is too slow!
    }

    for ( auto p: primitives ) {
        worldBound = Union(worldBound, p->BB());
    }
    for ( auto light: lights )
        light->Preprocess(* this);

}
