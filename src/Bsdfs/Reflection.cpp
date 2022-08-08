#include "Reflection.hpp"
#include "spdlog/spdlog.h"
#include "../Common/Frame.hpp"

Spectrum Bsdf::f(const vec3 & wo, const vec3 & wi, BxDFType flags) const {
    if (wo.z == 0) return Spectrum();
    bool reflect = Frame::cosTheta(wo) * Frame::cosTheta(wi)> 0;
    Spectrum f(0.f);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags) &&
            ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
             (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
            f += bxdfs[i]->f(wo, wi);
    return f;
}


int Bsdf::NumComponents(BxDFType flags) const {
    return 0;
}


Spectrum Lambertain::f(const vec3 & wo,const vec3 & wi) const {

//    if(wo.z<0 || wi.z<0){
//        return Spectrum();
//    }
    return albedo * Constant::INV_PI;
}

Lambertain::Lambertain(Spectrum & albedo) : albedo(albedo){

}

void Lambertain::LogInfo( ) const {
    spdlog::info("{0} albedo {1} {2} {3}",name,albedo.x,albedo.y,albedo.z);
}


