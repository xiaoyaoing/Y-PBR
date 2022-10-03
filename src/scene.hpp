#pragma  once

#include "nlohmann/json.hpp"
#include "Primitives/Primitive.hpp"
#include "Lights/Light.hpp"
#include "Hierancy/BVH.hpp"
#include "Primitives/EmbreeUtils.hpp"

class Light;
class Primitive;

class Scene {
public:
    Scene(const nlohmann::json j);

    // Scene and Light Intersection - Records the most recent intersection
    std::optional<Intersection>  intersect(const Ray& ray) const;

    //Determines whether the light and scene intersect - tests for occlusion
    bool intersectP(const Ray & ray) const ;
    //log some info for debug
    void logDebugInfo();
    void setUp();
    Bounds3 getWorldBound() const {
        return worldBound;
    }

    std::shared_ptr<BSDF> fetchBSDF(const std::string & bsdfName) const {
        if(bsdfs.contains(bsdfName)){
            return bsdfs.at(bsdfName);
        }
        return  bsdfs.at("default");
    }

    std::shared_ptr<BSDF> fetchBSDFFromJson(const nlohmann::json & bsdfJson){
        if(bsdfJson.is_string()){
            return fetchBSDF(bsdfJson.get <std::string>());
        }
        if(bsdfJson.is_object()){
            return fetchBSDF("default");
        }
    }

public :
    std::vector<std::shared_ptr<Light>> lights;
protected:
    /// add light if exists
    /// \param p the   object json
    /// \param l left  bound  idx of the
    /// \param r right bound  idx of the object
    void handleAddLight(const nlohmann::json & j,int l,int r);

    std::vector<std::shared_ptr<Primitive>> primitives;
    std::unordered_map<std::string,std::shared_ptr<BSDF>>  bsdfs;

    std::unique_ptr<BVHAccel> bvh;
    Bounds3 worldBound;
    RTCScene  _scene = nullptr;
    bool _useBVH;


    /*** debug variables ***/
};




