#include "Light.hpp"
#include "Texture/BitMapTexture.hpp"
#include "Common/Json.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"

class InfinteSphere : public  Infinite {
public:
    InfinteSphere(const std::shared_ptr<BitMapTexture<Spectrum>> emssision,const mat4 & toWorld) :
                _emission(emssision),_toWorld(toWorld),_toLocal(glm::transpose(toWorld)){}

    InfinteSphere(const Json & json);
    Spectrum
    sampleLi(const vec3 & ref, const vec2 & u, vec3 * wi, Float * pdf, Float * distance) const override;
    PositionAndDirectionSample sampleDirect(const vec2 & positionSample, const vec2 & u2) const override;
    Spectrum Le(const Ray & ray) const override;


    vec2 directionToUV(const vec3 &wi) const;
    vec2 directionToUV(const vec3 &wi, float &sinTheta) const;
    vec3 uvToDirection(vec2 uv, float &sinTheta) const;

    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;

    void logDebugInfo() const ;
private:
    Spectrum Power( ) override;
    void Preprocess(const Scene & scene) override;
protected:
    std::shared_ptr<BitMapTexture<Spectrum>> _emission;
    mat4 _toWorld;
    mat4 _toLocal;
};

