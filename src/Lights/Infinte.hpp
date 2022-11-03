#include "Light.hpp"
#include "Texture/BitMapTexture.hpp"



class InfinteSphere : public  Light {
public:
    InfinteSphere(const std::shared_ptr<BitMapTexture<Spectrum>> emssision,const mat4 & toWorld) :
                 Light((int)LightFlags::Infinite),_emission(emssision),_toWorld(toWorld),_toLocal(glm::transpose(toWorld)){}

    Spectrum
    sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;
    LightSampleResult sampleDirect(const vec2 & positionSample, const vec2 & u2) override;
    Spectrum directLighting(const Intersection & intr) const override;
    Spectrum environmentLighting(const Ray & ray) const override;


    vec2 directionToUV(const vec3 &wi) const;
    vec2 directionToUV(const vec3 &wi, float &sinTheta) const;
    vec3 uvToDirection(vec2 uv, float &sinTheta) const;

    Float directPdf(const Intersection & pShape, const vec3 & ref) const override;

    void logDebugInfo() const ;
private:
    Spectrum Power( ) override;
    void Preprocess(const Scene & scene) override;
protected:
    std::shared_ptr<BitMapTexture<Spectrum>> _emission;
    vec3 _worldCenter;
    Float _worldRadius;
    mat4 _toWorld;
    mat4 _toLocal;
};

