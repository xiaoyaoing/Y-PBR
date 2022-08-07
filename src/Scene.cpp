
#include "scene.hpp"
#include "iostream"
#include "Bsdfs/Bsdf.hpp"
#include "Bsdfs/BsdfFactory.hpp"
#include "spdlog/spdlog.h"

#include "Primitives/Triangle.hpp"
#include "Primitives/Sphere.hpp"
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
    auto loadQuadric = [](const nlohmann::json & json,std::shared_ptr<Bsdf> bsdf){
        return nullptr;
    };

    std::unordered_map<std::string,
        std::function<std::shared_ptr<Primitive>(const nlohmann::json,std::shared_ptr<Bsdf> bsdf)>
        >loadMap{{"object",loadObject},
                    {"sphere",loadSphere},
                    {"triangle",loadTriangle},
                    {"quadric",loadQuadric}};


    std::unordered_map<std::string,
                       std::shared_ptr<TriangleMesh>> mesh_set;  //mesh-set

    auto vertices =
            getOptional(j, "vertices", std::unordered_map<std::string, std::vector<vec3>>());  //vertices_map

    for(auto p : j.at("primitives")){
        //get bsdf
        std::string bsdf_str=p.at("bsdf");
        std::shared_ptr<Bsdf> bsdf = bsdfs[bsdf_str];
        if(!bsdf.get()){
            bsdf = bsdfs["default"];
        }

        //get transform
        std::unique_ptr<Transform> transform;
        if (p.find("position") != p.end() || p.find("scale") != p.end() || p.find("rotation") != p.end())
        {
            transform = std::make_unique<Transform>(
                    getOptional(p, "position", vec3(0.0)),
                    getOptional(p, "scale", vec3(1.0)),
                    radians(getOptional(p, "rotation", vec3(0.0)))
            );
        }


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
                            CreateTriangleMesh(transform.get(),&v,nullptr,nullptr,nullptr,bsdf);

                }
                std::shared_ptr<TriangleMesh> mesh=mesh_set[meshName];
                triangles_v = p.at("triangles").get<std::vector<std::vector<size_t>>>();
                triangles   = getTrianglesFromMesh(mesh,triangles_v);
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

        if(primitive)
            if(transform) primitive->transform(*transform);
            primitives.push_back(primitive);
        handleAddLight(p,primitives.size()-1,primitives.size());
        }
    }

    spdlog::info("{} Primitives",primitives.size());
    spdlog::info("{} lights",lights.size());
    spdlog::info("{} Bsdfs",bsdfs.size());
}


void Scene::handleAddLight(const nlohmann::json & p,size_t l,size_t r){
    if(p.find("light")!=p.end()){
        auto lightJson=p.at("light");
        std::string lightType=lightJson.at("type");

        if(lightType=="area"){
            Spectrum  albedo = lightJson.at("albedo");
            for(size_t i=l;i<r;i++){
                auto light = std::make_shared <AreaLight>(this->primitives[i],albedo);
                this->primitives[i]->areaLight =light;
                lights.push_back(light);
            }
        }
    }
    return ;
}



std::optional<Intersection> Scene::intersect(const Ray &ray) const {
    std::optional<Intersection> minIntersection;
    Ray _ray(ray);
    for(auto primitive : primitives) {
//        if(primitive->bsdf.get()== nullptr){
//            int k=1;
//        }
        auto its = primitive->intersect(_ray);
        if ( its.has_value() )
        {
            minIntersection = its;
        }
    }



//    if(minIntersection.has_value()){
//
////        minIntersection->n = minIntersection->primitive->n
//    }
//
    return minIntersection;
}
