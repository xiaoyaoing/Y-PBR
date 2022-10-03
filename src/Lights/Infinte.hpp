#include "Light.hpp"
#include "Texture/BitMapTexture.hpp"



class InfinteSphere : public  Light {
public:
    InfinteSphere(const std::shared_ptr<BitMapTexture> emssision) :
                 Light((int)LightFlags::Infinite),_emission(emssision){}
    Spectrum
    Sample_Li(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;
    Spectrum directLighting(const Intersection & intr) const override;

    vec2 directionToUV(const vec3 &wi) const;
    vec2 directionToUV(const vec3 &wi, float &sinTheta) const;
    vec3 uvToDirection(vec2 uv, float &sinTheta) const;

    Float directPdf(const Intersection & pShape, const vec3 & ref) const override;
    void logDebugInfo() const ;

    Spectrum environmentLighting(const Ray & ray) const override;

private:
    Float Power( ) override;

    void Preprocess(const Scene & scene) override;

protected:
    std::shared_ptr<BitMapTexture> _emission;
    vec3 _worldCenter;
    Float _worldRadius;
};

