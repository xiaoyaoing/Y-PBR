#include "Hair.h"

const int pMax = 3;

static Float I0(Float x){

}

Spectrum Hair::f(const SurfaceScatterEvent & event) const {
    Float sinThetaO = event.wo.y;
    Float cosThetaO = sqrt(1-sinThetaO * sinThetaO);
    Float phiO = std::atan2(event.wo.x,event.wo.z);
    Float thetaO = std::asin(clamp(sinThetaO,-1.0,1.0));

    Float sinThetaI = event.wi.x;
    Float cosThetaI = sqrt(1-sinThetaI * sinThetaI);
    Float thetaI = std::asin(clamp(sinThetaI,-1.0,1.0));

    Float thetaD = (thetaO-thetaI) * 0.5;
    Float cosThetaD = cos(thetaD);

    Float thetaOR = thetaO - 2.0 * _scaleAngel;
    Float thetaOTT = thetaO + _scaleAngel;
    Float thetaOTRT = thetaO + 4 * _scaleAngel;

    Float MR = M(_vR, sin(thetaOR), sinThetaI, cos(thetaOR), cosThetaI);
    Float MTT = M(_vR, sin(thetaOTT), sinThetaI, cos(thetaOTT), cosThetaI);
    Float MTRT = M(_vR, sin(thetaOTRT), sinThetaI, cos(thetaOTRT), cosThetaI);


    return Spectrum();
}

Float Hair::Pdf(const SurfaceScatterEvent & event) const {
    return 0;
}

Spectrum Hair::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    return Spectrum();
}

void Hair::LogInfo( ) const {

}

Float Hair::M(Float v, Float sinThetaO, Float sinThetaI, Float cosThetaO, Float cosThetaI) const {
    Float a = -sinThetaO * sinThetaI / v;
    Float b = cosThetaI * cosThetaO / v;
    Float csch = 2/(exp(1/v)-exp(-1/v));

    return csch * exp(a) * I0(b);
}
