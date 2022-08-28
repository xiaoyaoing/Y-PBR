
#include "scene.hpp"
#include "iostream"
#include "Bsdfs/Reflection.hpp"
#include "Bsdfs/BsdfFactory.hpp"
#include <spdlog/spdlog.h>

#include "Primitives/Triangle.hpp"
#include "Primitives/Sphere.hpp"
#include "Primitives/Quad.hpp"
#include "Primitives/Cube.hpp"

static int occludedCount =0 ;
static int occludedCount1 =0;


Scene::Scene(const nlohmann::json j) {
//    std::unordered_map<std::string,nlohmann::json> bsdf_jsons =j.at("bsdfs");
//
//    std::vector<nlohmann::json> bsdf_jsons = j.at("bsdfs");
//
//    std::unordered_map<std::string,std::shared_ptr<Bsdf>> bsdfs;
    auto bsdfs = BsdfFactory::LoadBsdfsFromJson(j.at("bsdfs"));
//    for(auto item:bsdf_jsons){
//        bsdfs[item.first]=BsdfFactory::LoadBsdfFromJson(item.second);
//    }
    auto loadObject = [](const nlohmann::json & json,std::shared_ptr<Bsdf> bsdf){
        return nullptr;
    };
    auto loadSphere = [](const nlohmann::json & json,std::shared_ptr<Bsdf> bsdf){
        return std::make_shared<Sphere>(json.at("radius"),bsdf);
    };
    auto loadTriangle = [](const nlohmann::json & json,std::shared_ptr<Bsdf> bsdf){
        return nullptr;
    };
    auto loadQuad = [](const nlohmann::json & json,std::shared_ptr<Bsdf> bsdf){
        return std::make_shared<Quad>(bsdf);
    };
    auto loadCube = [](const nlohmann::json & j,std::shared_ptr<Bsdf> bsdf){
        return std::make_shared <Cube>(bsdf);
    };


    std::unordered_map<std::string,
        std::function<std::shared_ptr<Primitive>(const nlohmann::json,std::shared_ptr<Bsdf> bsdf)>
        >loadMap{{"object",loadObject},
                    {"sphere",loadSphere},
                    {"triangle",loadTriangle},
                    {"quad",loadQuad},
                    {"cube",loadCube}
    };


    std::unordered_map<std::string,
                       std::shared_ptr<TriangleMesh>> mesh_set;  //mesh-set

    auto vertices =
            getOptional(j, "vertices", std::unordered_map<std::string, std::vector<vec3>>());  //vertices_map

    for(auto p : j.at("primitives")){
        //get bsdf
        std::shared_ptr<Bsdf> bsdf;
        std::string bsdf_str=p.at("bsdf");
        if(!bsdfs.contains(bsdf_str)){
            spdlog::error("{0},bsdf not found",bsdf_str);
            bsdf = bsdfs["default"];
        }
        else {
            bsdf=bsdfs[bsdf_str];
        }

        //get transform
        std::unique_ptr<Transform> transform;
        if(p.contains("transform"))
            transform = std::make_unique<Transform>(p["transform"]);


        std::string type = p.at("type");

        if(type=="object"){  // load shapes
            std::vector<vec3> v, n;
            std::vector<std::vector<size_t>> triangles_v, triangles_vt, triangles_vn;
//            auto vertices = getOptional(j, "vertices", std::unordered_map<std::string,vec3 []>());
            std::vector<std::shared_ptr<Primitive>> triangles;

            if(p.find("vertex_set")!=p.end()){//vertex object
                std::string meshName =p.at("vertex_set");

                if(!mesh_set.contains(meshName)){
                    v = vertices.at(p.at("vertex_set"));
                    mesh_set[meshName]=
                            CreateTriangleMesh(transform.get(),&v,nullptr,nullptr,nullptr);

                }
                std::shared_ptr<TriangleMesh> mesh=mesh_set[meshName];
                triangles_v = p.at("triangles").get<std::vector<std::vector<size_t>>>();
                triangles   = getTrianglesFromMesh(mesh,triangles_v,bsdf);
                }
            else{                                //.obj
                continue;
               }
            primitives.insert(primitives.end(),triangles.begin(),triangles.end());
            handleAddLight(p,primitives.size()-triangles.size(),primitives.size());

        }
        else  // load single shape
        {
        std::shared_ptr<Primitive>  primitive= loadMap[type](p,bsdf);

        if(transform) primitive->transform(*transform);
        if(//primitive->bsdf->name=="light" || primitive->bsdf->name=="shortBox" || primitive->bsdf->name=="tallBox"
         type!="sphere" || true
      // type=="quad"
        )
        {    primitives.push_back(primitive);
             handleAddLight(p,primitives.size()-1,primitives.size());
        }
        }
    }

//    primitives=std::vector<std::shared_ptr<Primitive>>(primitives.begin(),primitives.begin()+2);
    for(auto i:primitives){
        i->bsdf->LogInfo();
    }

    spdlog::info("{} Primitives",primitives.size());
    spdlog::info("{} lights",lights.size());
    spdlog::info("{} Bsdfs",bsdfs.size());

    getOptional(j["renderer"],"scene_bvh",_useBVH);

    if(_useBVH){
        bvh = std::make_unique<BVHAccel>(primitives);
    }


}


void Scene::handleAddLight(const nlohmann::json & p,size_t l,size_t r){
    if(p.find("light")!=p.end()){
        auto lightJson=p.at("light");
        std::string lightType=lightJson.at("type");

        if(lightType=="area"){
            Spectrum  albedo = lightJson.at("albedo");
            Float totalArea = 0,invArea;
            for(size_t i=l;i<r;i++) totalArea+=primitives[i]->Area();
            invArea=1/ totalArea;
            for(size_t i=l;i<r;i++){
                // flux to radiosity
                auto light = std::make_shared <AreaLight>(this->primitives[i],albedo * invArea) ;

                this->primitives[i]->areaLight =light;
                lights.push_back(light);
            }
        }
    }
    Spectrum  emission;
    if( containsAndGet(p,"emission",emission)){
        for(size_t i=l;i<r;i++){
            auto light = std::make_shared <AreaLight>(this->primitives[i],emission) ;
            this->primitives[i]->areaLight =light;
            lights.push_back(light);
        }
    }

    return ;
}



std::optional<Intersection> Scene::intersect(const Ray &ray) const {
    if(_useBVH){
        return bvh->intersect(ray);
    }

    std::optional<Intersection> minIntersection;
    Ray _ray(ray);
    for(auto primitive : primitives) {
//        if(primitive->bsdf.get()== nullptr){
//            int k=1;
//        }
        if(primitive->bsdf->name=="ceiling"){
            int k=1;
        }
        auto its = primitive->intersect(_ray);
        if ( its.has_value() )
        {
            minIntersection = its;
        }
    }
    return minIntersection;
}

bool Scene::intersectP(const Ray & ray) const {
    if(_useBVH ){
        //return bvh->intersectP(ray);
    }
    else
    {
        Ray _ray(ray);
        assert(_ray.farT ==ray.farT);
        for(auto primitive:primitives){
         auto its = primitive->intersect(_ray);
         if(its.has_value())
         {
             return true;
         }
         }
    }

    return false;
}

void Scene::logDebugInfo() {
    spdlog::info("Occluded {0} {1}",occludedCount,occludedCount1);
}
