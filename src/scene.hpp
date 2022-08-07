#pragma  once

#include "nlohmann/json.hpp"
#include "Primitives/Primitive.hpp"
#include "Lights/Light.hpp"

class Light;
class Primitive;

class Scene {
public:
    Scene(const nlohmann::json j);

    std::optional<Intersection>  intersect(const Ray& ray) const;

    std::vector<std::shared_ptr<Light>> lights;
private:
    /// add light if exists
    /// \param p the   object json
    /// \param l left  bound  idx of the
    /// \param r right bound  idx of the object
    void handleAddLight(const nlohmann::json & j,size_t l,size_t r);

    std::vector<std::shared_ptr<Primitive>> primitives;
    bool _useBVH;


//    std::unordered_map<std::string,std::shared_ptr<Bsdf>>  bsdfs;
};


