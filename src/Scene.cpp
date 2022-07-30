
#include "scene.hpp"
#include "iostream"
#include "Bsdfs/Bsdf.hpp"
Scene::Scene(const nlohmann::json j) {
    std::unordered_map<std::string,std::shared_ptr<Bsdf>> bsdfs = j.at("bsdfs");



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

    for(auto p : j.at("primitives")){
        //get bsdf
        std::string bsdf_str=p.at("bsdf");
        std::shared_ptr<Bsdf> bsdf = bsdfs[bsdf_str];

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

        }
        else  // load single shape
        {
        std::shared_ptr<Primitive>  primitive= loadMap[type](p,bsdf);

        if(primitive)
            if(transform) primitive->transform(*transform);
            primitives.push_back(primitive);
        }
    }
}



bool Scene::intersect(const Ray &ray,Intersection & intersection) const {
    bool intersected;
    Ray _ray(ray);
    for(auto primitive : primitives){
        if(primitive->intersect(_ray,intersection))
            intersected = true;
    }
    return intersected;
}
