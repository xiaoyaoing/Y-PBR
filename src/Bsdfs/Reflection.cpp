#include "Reflection.hpp"
#include <spdlog/spdlog.h>
#include "../Common/Frame.hpp"
#include "../Sampler/Warp.hpp"
#include "Fresnel.hpp"

#include "ComplexIor.hpp"
//Spectrum Bsdf::f(const vec3 & wo, const vec3 & wi, BXDFType flags) const {
//    if ( wo.z == 0 ) return Spectrum();
//    bool reflect = Frame::cosTheta(wo) * Frame::cosTheta(wi) > 0;
//    Spectrum f(0.f);
//    for ( int i = 0 ; i < nBXDFs ; ++ i )
//        if ( BXDFs[i]->MatchesFlags(flags) &&
//             ( ( reflect && ( BXDFs[i]->m_type & BSDF_REFLECTION ) ) ||
//               ( ! reflect && ( BXDFs[i]->m_type & BSDF_TRANSMISSION ) ) ) )
//            f += BXDFs[i]->f(wo, wi);
//    return f;
//}
//
//
//int Bsdf::NumComponents(BXDFType flags) const {
//    int MatchCount = 0;
//    for ( int i = 0 ; i < nBXDFs ; i ++ ) {
//        if ( BXDFs[i]->MatchesFlags(flags) )
//            MatchCount ++;
//    }
//    return MatchCount;
//}
//
//
//Spectrum Bsdf::sampleF(const vec3 & wo, vec3 * wi, const vec2 & u, Float * pdf,
//                       BXDFType type,
//                       BXDFType * sampledType) const {
//    int matchingComps = NumComponents(type);
//
//    if ( matchingComps == 0 ) {
//        * pdf = 0;
//        if ( sampledType ) * sampledType = BXDFType(0);
//        return Spectrum(0);
//    }
//
//    int comp =
//            std::min((int) std::floor(u[0] * matchingComps), matchingComps - 1);
//
//    BXDF * bxdf = nullptr;
//    int count = comp;
//    for ( int i = 0 ; i < nBXDFs ; ++ i )
//        if ( BXDFs[i]->MatchesFlags(type) && count -- == 0 ) {
//            bxdf = BXDFs[i];
//            break;
//        }
//
//    Spectrum f = bxdf->sampleF(wo, wi, u, pdf, sampledType);
//    if ( * pdf == 0 ) {
//        if ( sampledType ) * sampledType = BXDFType(0);
//        return Spectrum();
//    }
//
//    // avoid calculate Specular pdf directly
//    if ( ! ( bxdf->m_type & BSDF_SPECULAR ) && matchingComps > 1 )
//        for ( int i = 0 ; i < nBXDFs ; ++ i )
//            if ( BXDFs[i] != bxdf && BXDFs[i]->MatchesFlags(type) )
//                * pdf += BXDFs[i]->Pdf(wo, * wi);
//    if ( matchingComps > 1 ) * pdf /= matchingComps;
//
//    if ( ! ( bxdf->m_type & BSDF_SPECULAR ) ) {
//        bool reflect = Frame::cosTheta(wo) * Frame::cosTheta(* wi) > 0;
//        f = Spectrum(0);
//        for ( int i = 0 ; i < nBXDFs ; ++ i )
//            if ( BXDFs[i]->MatchesFlags(type) &&
//                 ( ( reflect && ( BXDFs[i]->m_type & BSDF_REFLECTION ) ) ||
//                   ( ! reflect && ( BXDFs[i]->m_type & BSDF_TRANSMISSION ) ) ) )
//                f += BXDFs[i]->f(wo, * wi);
//    }
//
//    return f;
//}


Spectrum LambertainR::f(const SurfaceScatterEvent & event) const {

    if(event.wo.z<0|| event.wi.z<0){
        return Spectrum();
    }
    Spectrum  albedo = m_albedo->Evaluate(event.its);
    return albedo * Constant::INV_PI ;
}



void LambertainR::LogInfo( ) const {
    Spectrum  albedo = m_albedo->Evaluate();
    spdlog::info("{0} albedo {1} {2} {3}", "LambertainReflection", albedo.x, albedo.y, albedo.z);
}


Spectrum LambertainR::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    
     event.wi=Warp::squareToUniformHemisphere(u);
     event.pdf=Warp::squareToUniformHemispherePdf(event.wi);
     event.sampleType=BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);

     return m_albedo->Evaluate(event.its) * Constant::INV_PI;
}


Spectrum LambertainT::f(const SurfaceScatterEvent & event) const {
    return m_albedo->Evaluate(event.its) * Constant::INV_PI;
}

void LambertainT::LogInfo( ) const {

}



Spectrum LambertainT::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {
    return Spectrum();
}




Float LambertainR::Pdf(const SurfaceScatterEvent & event) const {
    if (!MatchesFlags(event.sampleType))
        return 0.0f;
    if (event.wi.z <= 0.0f || event.wo.z <= 0.0f)
        return 0.0f;
    return Warp::squareToUniformHemispherePdf(event.wo);
}

Spectrum SpecularR::f(const SurfaceScatterEvent & event) const {
    const vec3 & wi = event.wi;
    const vec3 & wo = event.wo;
    if(Frame::Reflect(wi)==wo)
        return Spectrum(1.0);

    return Spectrum(0.0);
}

Float SpecularR::Pdf(const SurfaceScatterEvent & event) const {
    return 0;
}


void SpecularR::LogInfo( ) const {
    spdlog::info("Specular Reflection");
}



Spectrum
SpecularR::sampleF(SurfaceScatterEvent & event, const vec2 & u) const {

    event.wi= Frame::Reflect(event.wo);
    event.pdf=1;
    event.sampleType = BXDFType(BSDF_SPECULAR | BSDF_REFLECTION);
    return Spectrum(1.0) / abs(event.wi.z);
}

