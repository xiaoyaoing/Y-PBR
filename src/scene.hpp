#pragma  once
//2022/7/13

#include "nlohmann/json.hpp"
#include "Primitives/Primitive.hpp"

class Scene {
public:
    Scene(const nlohmann::json j);

    bool  intersect(const Ray& ray,Intersection & intersection) const;
private:
    std::vector<std::shared_ptr<Primitive>> primitives;
//    std::unordered_map<std::string,std::shared_ptr<Bsdf>>  bsdfs;
};


