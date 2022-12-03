#pragma  once

#include "Common/Json.hpp"
#include "Primitives/Primitive.hpp"
#include "Lights/Light.hpp"
#include "Hierancy/BVH.hpp"
#include "Primitives/EmbreeUtils.hpp"

class Light;
class Primitive;
class Medium;

struct RenderOptions{
    int spp;
    int sppStep;
    int maxBounces;
    RenderOptions(const Json & renderJson){
        spp = getOptional(renderJson,"spp",64);
        sppStep = getOptional(renderJson,"spp_step",1);
    }
};

class Scene {
public:
    Scene(const Json sceneJson);

    // Scene and Light Intersection - Records the most recent intersection
    std::optional<Intersection>  intersect(Ray & ray) const;

    //Determines whether the light and scene intersect - tests for occlusion
    bool intersectP(const Ray & ray) const ;
    //log some info for debug
    void logDebugInfo();
    void build();
    Bounds3 getWorldBound() const {
        return worldBound;
    }

    std::shared_ptr<BSDF> fetchBSDF(const std::string & bsdfName) const;
    std::shared_ptr<Medium> fetchMedium(const std::string & mediumName) const;

    void AddPrimitive(std::shared_ptr<Primitive> prim) {primitives.push_back(prim);}
    void AddLight(std::shared_ptr<Light> light) {lights.push_back(light);}
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
    std::unordered_map<std::string,std::shared_ptr<Medium>>  mediums;


    Bounds3 worldBound;
    RTCScene  _scene = nullptr;

    bool _useBVH;
    std::unique_ptr<BVHAccel> bvh;
};




