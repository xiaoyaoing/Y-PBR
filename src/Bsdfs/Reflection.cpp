#include "Reflection.hpp"
#include <spdlog/spdlog.h>
#include "../Common/Frame.hpp"
#include "../Sampler/Warp.hpp"
#include "Fresnel.hpp"
Spectrum Bsdf::f(const vec3 & wo, const vec3 & wi, BXDFType flags) const {
    if ( wo.z == 0 ) return Spectrum();
    bool reflect = Frame::cosTheta(wo) * Frame::cosTheta(wi) > 0;
    Spectrum f(0.f);
    for ( int i = 0 ; i < nBXDFs ; ++ i )
        if ( BXDFs[i]->MatchesFlags(flags) &&
             ( ( reflect && ( BXDFs[i]->m_type & BSDF_REFLECTION ) ) ||
               ( ! reflect && ( BXDFs[i]->m_type & BSDF_TRANSMISSION ) ) ) )
            f += BXDFs[i]->f(wo, wi);
    return f;
}


int Bsdf::NumComponents(BXDFType flags) const {
    int MatchCount = 0;
    for ( int i = 0 ; i < nBXDFs ; i ++ ) {
        if ( BXDFs[i]->MatchesFlags(flags) )
            MatchCount ++;
    }
    return MatchCount;
}


Spectrum Bsdf::sampleF(const vec3 & wo, vec3 * wi, const vec2 & u, Float * pdf,
                       BXDFType type,
                       BXDFType * sampledType) const {
    int matchingComps = NumComponents(type);

    if ( matchingComps == 0 ) {
        * pdf = 0;
        if ( sampledType ) * sampledType = BXDFType(0);
        return Spectrum(0);
    }

    int comp =
            std::min((int) std::floor(u[0] * matchingComps), matchingComps - 1);

    BXDF * bxdf = nullptr;
    int count = comp;
    for ( int i = 0 ; i < nBXDFs ; ++ i )
        if ( BXDFs[i]->MatchesFlags(type) && count -- == 0 ) {
            bxdf = BXDFs[i];
            break;
        }

    Spectrum f = bxdf->sampleF(wo, wi, u, pdf, sampledType);
    if ( * pdf == 0 ) {
        if ( sampledType ) * sampledType = BXDFType(0);
        return Spectrum();
    }

    // avoid calculate Specular pdf directly
    if ( ! ( bxdf->m_type & BSDF_SPECULAR ) && matchingComps > 1 )
        for ( int i = 0 ; i < nBXDFs ; ++ i )
            if ( BXDFs[i] != bxdf && BXDFs[i]->MatchesFlags(type) )
                * pdf += BXDFs[i]->Pdf(wo, * wi);
    if ( matchingComps > 1 ) * pdf /= matchingComps;

    if ( ! ( bxdf->m_type & BSDF_SPECULAR ) ) {
        bool reflect = Frame::cosTheta(wo) * Frame::cosTheta(* wi) > 0;
        f = Spectrum(0);
        for ( int i = 0 ; i < nBXDFs ; ++ i )
            if ( BXDFs[i]->MatchesFlags(type) &&
                 ( ( reflect && ( BXDFs[i]->m_type & BSDF_REFLECTION ) ) ||
                   ( ! reflect && ( BXDFs[i]->m_type & BSDF_TRANSMISSION ) ) ) )
                f += BXDFs[i]->f(wo, * wi);
    }

    return f;
}


Spectrum LambertainR::f(const vec3 & wo, const vec3 & wi) const {

//    if(wo.z<0 || wi.z<0){
//        return Spectrum();
//    }
    return albedo * Constant::INV_PI
         ;
}

LambertainR::LambertainR(Spectrum & albedo) : albedo(albedo), useCosineSample(true),
                                              BXDF(BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE))

{

}

void LambertainR::LogInfo( ) const {
    spdlog::info("{0} albedo {1} {2} {3}", "LambertainReflection", albedo.x, albedo.y, albedo.z);
}


Spectrum LambertainR::sampleF(const vec3 & wo, vec3 * wi,
                              const vec2 & u, Float * pdf,
                              BXDFType * sampledType) const {
    if(useCosineSample){
        *wi=Warp::squareToCosineHemisphere(u);
        *pdf=Warp::squareToCosineHemispherePdf(*wi);
    }
    else{
        *wi=Warp::squareToUniformHemisphere(u);
        *pdf=Warp::squareToUniformHemispherePdf(*wi);
    }

    //*wi = vec3(0,0,1);

    *sampledType=BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);

    return f(wo,*wi);
}


Spectrum LambertainT::f(const vec3 & wo, const vec3 & wi) const {
    return albedo * Constant::INV_PI;
}

void LambertainT::LogInfo( ) const {

}

LambertainT::LambertainT(Spectrum & albedo) : albedo(albedo),
                                              BXDF(BXDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE)) {

}

Spectrum LambertainT::sampleF(const vec3 & wo, vec3 * wi, const vec2 & u, Float * pdf, BXDFType * sampledType) const {
    return Spectrum();
}

void Bsdf::LogInfo( ) {

}


Float LambertainR::Pdf(const vec3 & wo, const vec3 & wi) const {
    return 0;
}

Spectrum SpecularR::f(const vec3 & wo, const vec3 & wi) const {
    if(Frame::Reflect(wi)==wo)
        return Spectrum(1.0);

    return Spectrum(0.0);
}

Float SpecularR::Pdf(const vec3 & wo, const vec3 & wi) const {
    return 0;
}


void SpecularR::LogInfo( ) const {
    spdlog::info("Specular Reflection");
}

SpecularR::SpecularR( ) : BXDF(BXDFType(BSDF_REFLECTION | BSDF_SPECULAR)) {

}

Spectrum
SpecularR::sampleF(const vec3 & wo, vec3 * wi,
                            const vec2 & u, Float * pdf,
                            BXDFType * sampledType) const {

    *wi= Frame::Reflect(wo);
    *pdf=1;
    *sampledType = BXDFType(BSDF_SPECULAR | BSDF_REFLECTION);
    return Spectrum(1.0);
}

Spectrum Dielectric::f(const vec3 & wo, const vec3 & wi) const {
   return Spectrum(0);
}

Float Dielectric::Pdf(const vec3 & wo, const vec3 & wi) const {
    // avoid compute specular pdf
    return 0;
}

Spectrum Dielectric::sampleF(const vec3 & wo, vec3 * wi, const vec2 & u, Float * pdf, BXDFType * sampledType) const {
    bool sampleT = enableT;
    Float  eta = wo.z < 0.0f ? ior : invIor;
    Float cosThetaT;
    Float F=Fresnel::DielectricReflectance(eta, std::abs(wo.z), cosThetaT);



    Float reflectionProbability=sampleT?F:1;

    if(u[0]<=reflectionProbability) {   //reflection
        *wi =Frame::Reflect(wo);
        *pdf = reflectionProbability;
        *sampledType = BXDFType(BSDF_REFLECTION | BSDF_SPECULAR);
        return albedo * F;
    }
    else {
        if(reflectionProbability == 1.0){
            * pdf =0;
            return  Spectrum(0);
        }

        *wi = vec3(-eta*wo.x,-eta*wo.y,-std::copysign(cosThetaT,wo.z));
       // *wi = vec3(0,0,-std::copysign(1,wo.z));
        *pdf = 1-reflectionProbability;
        *sampledType = BXDFType(BSDF_TRANSMISSION | BSDF_SPECULAR);
        Spectrum  f = albedo * (1-F);
        return f;
    }

    return Spectrum();
}

void Dielectric::LogInfo( ) const {
    spdlog::info("{DielectricBXDF ior:{0} albedo:{1} enableT{2}:",ior, toColorStr(albedo),enableT);
}

Spectrum Conductor::f(const vec3 & wo, const vec3 & wi) const {
   return Spectrum();
}

Float Conductor::Pdf(const vec3 & wo, const vec3 & wi) const {
    return 0;
}

Spectrum Conductor::sampleF(const vec3 & wo, vec3 * wi, const vec2 & u, Float * pdf, BXDFType * sampledType) const {
    * wi = Frame::Reflect(wo);
    *pdf = 1;
    Spectrum  f = m_albedo *  Fresnel::conductorReflectance(m_eta, m_k, wo.z);
    *sampledType = m_type;
    return f;
}

void Conductor::LogInfo( ) const {

}
