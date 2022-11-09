#include <spdlog/spdlog.h>
#include "iostream"

#include "scene.hpp"
#include "IO/MeshIO.hpp"
#include "IO/FileUtils.hpp"

#include "Bsdfs/Reflection.hpp"
#include "Bsdfs/BsdfFactory.hpp"

#include "Primitives/TriangleMesh.hpp"
#include "Primitives/Sphere.hpp"
#include "Primitives/Quad.hpp"
#include "Primitives/Cube.hpp"

#include "Lights/Infinte.hpp"

#include "Texture/TextureFactory.hpp"
#include "Common/Transform.hpp"

static int occludedCount = 0;
static int occludedCount1 = 0;


Scene::Scene(const Json sceneJson) : options(RenderOptions(sceneJson.at("renderer"))){
    bsdfs = BSDFFactory::LoadBsdfsFromJson(sceneJson.at("bsdfs"));


    auto loadSphere = [](const Json & json, std::shared_ptr < BSDF > bsdf) {
        Float radius = getOptional(json,"radius",1.f);
        return std::make_shared < Sphere >(radius, bsdf);
    };

    auto loadQuad = [](const Json & json, std::shared_ptr < BSDF > bsdf) {
        return std::make_shared < Quad >(bsdf);
    };
    auto loadCube = [](const Json & j, std::shared_ptr < BSDF > bsdf) {
        return std::make_shared < Cube >(bsdf);
    };
    auto loadInfiniteSphere = [](const Json & j, std::shared_ptr < BSDF > bsdf) {
        return nullptr;
    };


    std::unordered_map < std::string,
            std::function < std::shared_ptr < Primitive >(const Json, std::shared_ptr < BSDF > bsdf) >
    > loadMap{{"sphere",          loadSphere},
              {"quad",            loadQuad},
              {"cube",            loadCube},
              {"infinite_sphere", loadInfiniteSphere},
    };

    //load primiteives
    for ( auto p: sceneJson.at("primitives") ) {
        std::shared_ptr < BSDF > bsdf = fetchBSDFFromJson(getOptional(p,"bsdf",std::string("null")));

        mat4 transform = getOptional(p,"transform",getIndentifyTransform());


        std::string type = p.at("type");
        spdlog::info(type);
        if ( type == "mesh" ) {
            std::shared_ptr < TriangleMesh > mesh = std::make_shared < TriangleMesh >();
            mesh->Load(p, * this, transform);
            primitives.push_back(mesh);
            // primitives.insert(primitives.end(),triangles.begin(),triangles.end());
            handleAddLight(p, primitives.size() - 1, primitives.size());
        } else {
            if(!loadMap.contains(type)) continue;
            std::shared_ptr < Primitive > primitive = loadMap[type](p, bsdf);
            if ( primitive ) primitive->transform(transform);
            if ( primitive //&& type!="cube"
            ) {
                primitives.push_back(primitive);
                handleAddLight(p, primitives.size() - 1, primitives.size());
            } else {
                if ( type == "infinite_sphere" ) {
                    auto bitMap = std::make_shared < BitMapTexture<Spectrum> >
                            (FileUtils::WorkingDir + p.at("emission").get < std::string >());
                    bitMap->LoadResources();
                    mat4 toWorld = getOptional(p,"transform",getIndentifyTransform());
                    lights.push_back(std::make_shared <InfinteSphere>(bitMap,toWorld));
                }
            }
        }
    }

//    primitives=std::vector<std::shared_ptr<Primitive>>(primitives.begin(),primitives.begin()+2);

    spdlog::info("{} Primitives", primitives.size());
    spdlog::info("{} lights", lights.size());
    spdlog::info("{} Bsdfs", bsdfs.size());

    _useBVH = getOptional(sceneJson["renderer"], "scene_bvh", true);
   // _useBVH = false;

}


void Scene::handleAddLight(const Json & p,int l,int r) {
//    if ( p.find("light") != p.end() ) {
//        auto lightJson = p.at("light");
//        std::string lightType = lightJson.at("type");
//
//        if ( lightType == "area" ) {
//            Spectrum albedo = lightJson.at("albedo");
//            Float totalArea = 0, invArea;
//            for ( size_t i = l ; i < r ; i ++ ) totalArea += primitives[i]->Area();
//            invArea = 1 / totalArea;
//            for ( size_t i = l ; i < r ; i ++ ) {
//                // flux to radiosity
//                auto light = std::make_shared < AreaLight >(this->primitives[i], albedo * invArea);
//                this->primitives[i]->areaLight = light;
//                lights.push_back(light);
//            }
//        }
//    }
    if( contains(p,"emission")){
        std::shared_ptr<Texture<Spectrum>> emission =  TextureFactory::LoadTexture <Spectrum>(p,"emission",Spectrum(1));
        for ( size_t i = l ; i < r ; i ++ ) {
            auto light = std::make_shared < AreaLight >(this->primitives[i], emission);
            this->primitives[i]->areaLight = light;
            lights.push_back(light);
        }
    }

    else if( contains(p,"power")){
        std::shared_ptr<Texture<Spectrum>> power =  TextureFactory::LoadTexture <Spectrum>(p,"power",Spectrum(1));
        Float totalScale=0;
        for(size_t i = l ; i < r ; i ++){
            totalScale+=1 / primitives[i]->powerToRadianceScale();
        }
        totalScale = 1/totalScale;
        power->setScale(totalScale);

        for ( size_t i = l ; i < r ; i ++ ) {
            auto light = std::make_shared < AreaLight >(this->primitives[i], power);
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
        if ( rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID ) {
            return std::nullopt;
        }
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
        return shadowRay.tfar != ray.farT;
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


void Scene::build( ) {


    if ( _useBVH ) {
        _scene = rtcNewScene(EmbreeUtils::getDevice());
        for ( const auto & primitve: primitives ) {
            RTCGeometry geometry = primitve->initRTC();
            auto geom = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_USER);
            rtcEnableGeometry(geom);
            rtcSetGeometryUserPrimitiveCount(geom, 1);
            rtcSetGeometryUserData(geom, primitve.get());
            rtcSetGeometryBoundsFunction(geom, & EmbreeUtils::instanceBoundsFunc, nullptr);
            rtcSetGeometryIntersectFunction(geom, & EmbreeUtils::instanceIntersectFunc);
            rtcSetGeometryOccludedFunction(geom, & EmbreeUtils::instanceOccludedFunc);
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
