#include "Light.hpp"
#include "Texture/BitMapTexture.hpp"
#include "Common/Frame.hpp"


class InfinteSphereCap : public  Light {
public:
//    InfinteSphere(const std::shared_ptr<BitMapTexture<Spectrum>> emssision,const mat4 & toWorld) :
//            Light((int)LightFlags::Infinite),_emission(emssision),_toWorld(toWorld),_toLocal(glm::transpose(toWorld)){}
    InfinteSphereCap(const Json & json);

    Spectrum
    sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;
    LightSampleResult sampleDirect(const vec2 & positionSample, const vec2 & u2) override;
    Spectrum Le(const Ray & ray) const override;


    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;

    void logDebugInfo() const ;
private:
    Spectrum Power( ) override;
    void Preprocess(const Scene & scene) override;
protected:
    Spectrum _emission;

    vec3 _capDir;
    Float _capAngle;
    Float _cosCapAngle;

    vec3 _worldCenter;
    Float _worldRadius;

    Frame _capFrame;
};

