#include "Reflection.hpp"

class Hair : public  BSDF{
public:
    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;
protected:
    Float M(Float v,Float sinThetaO,Float sinThetaI,Float cosThetaO,Float cosThetaI) const ;

private:
    Float _scaleAngel;
    Float _roughness;
    Float _betaR,_betaTT,_betaTRT;
    Float _vR,_vTT,_vTRT;
};