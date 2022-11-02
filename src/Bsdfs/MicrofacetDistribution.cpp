#include "MicrofacetDistribution.hpp"
#include "Reflection.hpp"


const static std::unordered_map<std::string,MircroDistributionEnum> distribMap = {
        {"beckmann",Beckmann},
        {"ggx",GGx},
        {"trowbridge_reitz",TrowbridgeReitz}
};

void from_json(const Json & j, MircroDistributionEnum & type) {
    if(!j.is_string()){
        _ERROR("Should be a distribution string")
    }
    if(!distribMap.count(j))
        _ERROR("Not supported distribution");
    type =  distribMap.at(j);
}

MicrofacetDistribution::~MicrofacetDistribution( ) noexcept {}

Float MicrofacetDistribution::Pdf(const vec3 & wo, const vec3 & wh, const vec2 &alphaxy) const {
    if (sampleVisibleArea)
        return D(wh, alphaxy) * G1(wo,alphaxy) * absDot(wo, wh) / AbsCosTheta(wo);
    else
        return D(wh, alphaxy) * AbsCosTheta(wh);
}

Float BeckmannDistribution::roughnessToAlpha(float roughness) const {
    roughness = std::max(roughness, (Float)1e-3);
    return roughness;
    Float x = std::log(roughness);
    return 1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
}

Float BeckmannDistribution::D(const vec3 & wh, const vec2 & alphaxy) const {
    Float alphax =alphaxy.x;
    Float alphay =alphaxy.y;
    Float tan2Theta = Tan2Theta(wh);
    if (std::isinf(tan2Theta)) return 0.;
    Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
    return std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) +
                                  Sin2Phi(wh) / (alphay * alphay))) /
           (Constant::PI * alphax * alphay * cos4Theta);
}

Float BeckmannDistribution::Lambda(const vec3 & w, const vec2 & alphaxy) const {
    Float alphax =alphaxy.x;
    Float alphay =alphaxy.y;

    Float absTanTheta = std::abs(TanTheta(w));
    if (std::isinf(absTanTheta)) return 0.;
    // Compute _alpha_ for direction _w_
    Float alpha =
            std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
    Float a = 1 / (alpha * absTanTheta);
    if (a >= 1.6f) return 0;
    return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

vec3 BeckmannDistribution::Sample_wh(const vec3 & wo, const vec2 & u, const vec2 & alphaxy) const {
    Float alphax =alphaxy.x;
    Float alphay =alphaxy.y;

    //random sample  no importance sampling
    Float tan2Theta;
    Float phi;
    //todo add support sample visible area
    if(!sampleVisibleArea || true){
        if(alphax == alphay){
            tan2Theta = -alphax * alphax * std::log(1-u[0]);
            phi = u[1] * 2 * Constant::PI;
        }
        else {

        }
        Float cosTheta = std::sqrt(1 / (1+tan2Theta));
        Float sinTheta = std::sqrt(1 / (1+1/tan2Theta));
        vec3  wh = vec3(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);
        if(dot(wo,wh)<0) wh = -wh;
        return wh;
    }
    // see https://hal.inria.fr/hal-00996995v1/document // todo
    else {

    }
}

std::string BeckmannDistribution::ToString( ) const {
    return std::string();
}

Float TrowbridgeReitzDistribution::roughnessToAlpha(float roughness) const {
    return 0;
}

Float TrowbridgeReitzDistribution::D(const vec3 & wh, const vec2 & alphaxy) const {
    return 0;
}

Float TrowbridgeReitzDistribution::Lambda(const vec3 & w, const vec2 & alphaxy) const {
    return 0;
}

vec3 TrowbridgeReitzDistribution::Sample_wh(const vec3 & wo, const vec2 & u, const vec2 & alphaxy) const {
    return vec3();
}

std::string TrowbridgeReitzDistribution::ToString( ) const {
    return std::string();
}

