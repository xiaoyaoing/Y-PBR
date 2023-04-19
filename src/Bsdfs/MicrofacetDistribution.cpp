#include "MicrofacetDistribution.hpp"
#include "Reflection.hpp"


//const static std::unordered_map<std::string,MircroDistributionEnum> distribMap = {
//        {"beckmann",Beckmann},
//        {"ggx",GGx},
//        {"trowbridge_reitz",TrowbridgeReitz}
//};
//
//
//
//
//const static std::unordered_map<std::string,std::function<std::shared_ptr<MicrofacetDistribution>()>> _distribMap = {
//        {"beckmann",std::make_shared<Beckmann>},
//        {"ggx",std::make_shared<GGX>}
//};
//
//
//void from_json(const Json & j, MircroDistributionEnum & type) {
//    if(!j.is_string()){
//        _ERROR("Should be a distribution string")
//    }
//    if(!distribMap.count(j))
//        _ERROR("Not supported distribution");
//    type =  distribMap.at(j);
//}

MicrofacetDistribution::~MicrofacetDistribution( ) noexcept {}

Float MicrofacetDistribution::Pdf(const vec3 & wo, const vec3 & wh, const vec2 &alphaXY) const {
    if (sampleVisibleArea)
        return D(wh, alphaXY) * G1(wo,alphaXY) * absDot(wo, wh) / AbsCosTheta(wo);
    else
        return D(wh, alphaXY) * AbsCosTheta(wh);
}

Float Beckmann::roughnessToAlpha(Float roughness) const {
    roughness = std::max(roughness, (Float)1e-3);
    return roughness;
    Float x = std::log(roughness);
    return 1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
}

Float Beckmann::D(const vec3 & wh, const vec2 & alphaXY) const {
    Float alphax =alphaXY.x;
    Float alphay =alphaXY.y;
    Float tan2Theta = Tan2Theta(wh);
    if (std::isinf(tan2Theta)) return 0.;
    Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
    return std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) +
                                  Sin2Phi(wh) / (alphay * alphay))) /
           (Constant::PI * alphax * alphay * cos4Theta);
}

Float Beckmann::Lambda(const vec3 & w, const vec2 & alphaXY) const {
    Float alphax =alphaXY.x;
    Float alphay =alphaXY.y;

    Float absTanTheta = std::abs(TanTheta(w));
    if (std::isinf(absTanTheta)) return 0.;
    // Compute _alpha_ for direction _w_
    Float alpha =
            std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
    Float a = 1 / (alpha * absTanTheta);
    if (a >= 1.6f) return 0;
    return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

vec3 Beckmann::Sample_wh(const vec3 & wo, const vec2 & u, const vec2 & alphaXY) const {
    Float alphax =alphaXY.x;
    Float alphay =alphaXY.y;

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
            phi = std::atan(alphay / alphax *
                            std::tan(2 * Constant::PI * u[1] + 0.5 * Constant::PI));
            if (u[1] > 0.5) phi += Constant::PI;
            double sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            double alphaX2 = alphax * alphax, alphay2 = alphay * alphay;
            tan2Theta = -std::log(1-u[0]) /
                        (cosPhi * cosPhi / alphaX2 + sinPhi * sinPhi / alphay2);
        }
        Float cosTheta = std::sqrt(1 / (1+tan2Theta));
        Float sinTheta = std::sqrt(1 / (1+1/tan2Theta));
        vec3  wh = vec3(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);
        return wh;
    }
    // see https://hal.inria.fr/hal-00996995v1/document // todo
    else {
        _NOT_IMPLEMENT_ERROR
    }
}

std::string Beckmann::ToString( ) const {
    return std::string();
}

Float GGX::roughnessToAlpha(Float roughness) const {
    return roughness;

}

Float GGX::D(const vec3 & wh, const vec2 & alphaXY) const {
    Float alphaX = alphaXY.x;
    Float alphaY = alphaXY.y;
    Float ax2 = alphaX * alphaX;
    Float ay2 = alphaY * alphaY;
    vec3  wh2 = wh * wh;
    Float D = Constant::PI * alphaX *alphaY * pow(wh2.x/ax2+wh2.y/ay2+wh2.z,2);
    return 1/D;
}

Float GGX::Lambda(const vec3 & w, const vec2 & alphaXY) const {
    Float ax2 =  alphaXY.x * alphaXY.x;
    Float ay2 = alphaXY.y * alphaXY.y;
    vec3 v2 = w * w;
    Float Lambda = (-1 + sqrt(1 + (v2.x * ax2 + v2.y * ay2) / v2.z)) / 2;
    return Lambda;
}

vec3 GGX::Sample_wh(const vec3 & wo, const vec2 & u, const vec2 & alphaXY) const {
    Float alphaX = alphaXY.x,alphaY =alphaXY.y;
    if (wo.z < 0) {
        // Ensure the input is on top of the surface.
        return Sample_wh(-wo, u,alphaXY);
    }
    if(sampleVisibleArea){
        //see https://jcgt.org/published/0007/04/01/slides.pdf

        // Transform the incoming direction to the "hemisphere configuration".
        vec3 hemisphereDirOut= normalize(vec3(alphaX * wo.x, alphaY * wo.y, wo.z));
        // Parameterization of the projected area of a hemisphere.
        Float r = sqrt(u.x);
        Float phi = 2 * Constant::PI * u.y;
        Float t1 = r * cos(phi);
        Float t2 = r * sin(phi);
        // Vertically scale the position of a sample to account for the projection.
        Float s = (1 + hemisphereDirOut.z) / 2;
        t2 = (1 - s) * sqrt(1 - t1 * t1) + s * t2;
        // Point in the disk space
        vec3 diskN{t1, t2, sqrt(std::max(0.0f, 1 - t1*t1 - t2*t2))};
        // Reprojection onto hemisphere -- we get our sampled normal in hemisphere space.
//        vec3 T1 = normalize(vec3(-hemisphereDirOut.y,hemisphereDirOut.x,0));
//        vec3 T2 = cross(hemisphereDirOut,T1);
//        vec3 hemisphereN =  t1 * T1 + t2 * T2 + diskN.z * hemisphereDirOut;
        Frame frame(hemisphereDirOut);
        vec3 hemisphereN = frame.toWorld(vec3(t1,t2,diskN.z));
        // Transforming the normal back to the ellipsoid configuration
        auto wh =  normalize(vec3(alphaX * hemisphereN.x, alphaY *  hemisphereN.y, std::max(0.0f, hemisphereN.z)));
        if( hasNan(wh)){

        }
        return wh;
    }
    else {
        Float cosTheta, phi = (2 * Constant::PI) * u[1];
        if (alphaX == alphaY) {
            Float tanTheta2 = alphaX * alphaY * u[0] / (1.0f - u[0]);
            cosTheta = 1 / std::sqrt(1 + tanTheta2);
        } else {
            phi =
                    std::atan(alphaY / alphaX * std::tan(2 * Constant::PI * u[1] + .5f * Constant::PI));
            if (u[1] > .5f) phi += Constant::PI;
            Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            const Float alphaX2 = alphaX * alphaX, alphaY2 = alphaY * alphaY;
            const Float alpha2 =
                    1 / (cosPhi * cosPhi / alphaX2 + sinPhi * sinPhi / alphaY2);
            Float tanTheta2 = alpha2 * u[0] / (1 - u[0]);
            cosTheta = 1 / std::sqrt(1 + tanTheta2);
        }
        Float sinTheta =
                std::sqrt(std::max((Float )0., (Float )1. - cosTheta * cosTheta));
        vec3 wh = vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
        return  wh;
    }
}

std::string GGX::ToString( ) const {
    return std::string();
}

std::shared_ptr<MicrofacetDistribution>  LoadMicrofacetDistribution(const std::string & type){
    if(type == "beckmann")
        return std::make_shared <Beckmann>();
    if(type == "ggx")
        return std::make_shared <GGX>();
}


