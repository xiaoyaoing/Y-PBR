#include "Light.hpp"
#include "Texture/BitMapTexture.hpp"
#include "Common/Frame.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"


class InfinteSphereCap : public Infinite {
public:
//    InfinteSphere(const std::shared_ptr<BitMapTexture<Spectrum>> emssision,const mat4 & toWorld) :
//            Light((int)LightFlags::Infinite),_emission(emssision),_toWorld(toWorld),_toLocal(glm::transpose(toWorld)){}
    InfinteSphereCap(const Json & json);

    Spectrum
    sampleLi(const vec3 & ref, const vec2 & u, vec3 * wi, Float * pdf, Float * distance) const override;

    PositionAndDirectionSample sampleDirect(const vec2 & positionSample, const vec2 & dirSample) const override;

    Spectrum Le(const Ray & ray) const override;


    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;

    void logDebugInfo( ) const;

protected:
    Spectrum Power( ) override;

protected:
    Spectrum _emission;

    vec3 _capDir;
    Float _capAngle;
    Float _cosCapAngle;
    Frame _capFrame;
};

