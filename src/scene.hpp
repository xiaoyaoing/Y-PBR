#pragma  once

#include "Common/Json.hpp"
#include "Primitives/Primitive.hpp"
#include "Lights/Light.hpp"
#include "Hierancy/BVH.hpp"
#include "Primitives/EmbreeUtils.hpp"

class Light;
class Primitive;

struct RenderOptions{
    int spp;
    int maxBounces;
    RenderOptions(const Json & renderJson){
        spp = getOptional(renderJson,"spp",64);
    }
};

class Scene {
public:
    Scene(const Json sceneJson);

    // Scene and Light Intersection - Records the most recent intersection
    std::optional<Intersection>  intersect(const Ray& ray) const;

    //Determines whether the light and scene intersect - tests for occlusion
    bool intersectP(const Ray & ray) const ;
    //log some info for debug
    void logDebugInfo();
    void build();
    Bounds3 getWorldBound() const {
        return worldBound;
    }

    std::shared_ptr<BSDF> fetchBSDF(const std::string & bsdfName) const {
        if(bsdfs.contains(bsdfName)){
            return bsdfs.at(bsdfName);
        }
        return  bsdfs.at("default");
    }

    std::shared_ptr<BSDF> fetchBSDFFromJson(const Json & bsdfJson){
        if(bsdfJson.is_string()){
            return fetchBSDF(bsdfJson.get <std::string>());
        }
        if(bsdfJson.is_object()){
            return fetchBSDF("default");
        }
    }

    RenderOptions options;
public :
    std::vector<std::shared_ptr<Light>> lights;
protected:
    /// add light if exists
    /// \param p the   object json
    /// \param l left  bound  idx of the
    /// \param r right bound  idx of the object
    void handleAddLight(const Json & j,int l,int r);

    std::vector<std::shared_ptr<Primitive>> primitives;
    std::unordered_map<std::string,std::shared_ptr<BSDF>>  bsdfs;

    std::unique_ptr<BVHAccel> bvh;
    Bounds3 worldBound;
    RTCScene  _scene = nullptr;
    bool _useBVH;

    /*** debug variables ***/
};




