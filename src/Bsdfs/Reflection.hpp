#pragma  once

#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "Colors/Spectrum.hpp"
#include "BsdfTypes.hpp"
#include "Common/Texture.hpp"

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

//class BSDF;
////
//class Bsdf{
//
//public:
//    Spectrum f(const vec3 &wo , const vec3 &wi,
//               BXDFType flags = BSDF_ALL) const ;
//
//    ///sample a new dir
//    ///wi wo in local space
//    virtual Spectrum sampleF(const vec3 &wo, vec3 *wi, const vec2 &u,
//                             Float *pdf, BXDFType type,BXDFType *sampledType) const ;
//
//    int NumComponents(BXDFType flags = BSDF_ALL) const;
//
//    void Add(BXDF *b) {
//        BXDFs[nBXDFs++] = b;
//    }
//    void LogInfo();
//    std::string name;
//
////    Float eta(){
////        return BXDF[]
////    }
//private:
//    BXDF *BXDFs[8];
//    int nBXDFs = 0;
//    //for debug
//};

class BSDF{
public:

    BSDF(BXDFType type) : m_type(type),m_albedo(nullptr),m_bumpMap(nullptr) {};


    virtual Spectrum f(const SurfaceScatterEvent & event) const =0;

    virtual Float Pdf(const SurfaceScatterEvent & event) const =0;

    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const =0;


    virtual void LogInfo() const =0;

    bool MatchesFlags(BXDFType t) const {
        return ( m_type & t) !=0;
    }

    virtual  Float eta() const { return 1; }

    void setAlbedo(const std::shared_ptr<Texture<Spectrum>>  albedo){
        if(albedo == nullptr){
            ;
        }
        m_albedo = std::move(albedo);
    }


    void setBumpMap(const std::shared_ptr<Texture<Float>> &  bumpMap)  { m_bumpMap = bumpMap; }

    std::string name;
    std::shared_ptr<Texture<Spectrum>> m_albedo = nullptr ;
    int  temp;
protected:
    std::shared_ptr<Texture<Float>>    m_bumpMap = nullptr;
    BXDFType m_type;
};


class Mirror : public BSDF{

};

class LambertainR : public  BSDF{

public:
    LambertainR() : BSDF(BXDFType(BSDF_DIFFUSE | BSDF_REFLECTION)) { }

    virtual Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    virtual void LogInfo() const;


    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;
private:
    bool useCosineSample;
};

class LambertainT : public  BSDF{
public:

    virtual Spectrum f(const SurfaceScatterEvent & event) const override;

    virtual void LogInfo() const override;

    LambertainT() : BSDF(BXDFType(BSDF_DIFFUSE | BSDF_TRANSMISSION)) { }


    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

private:
};

class SpecularR : public BSDF{
public:


    virtual Spectrum f(const SurfaceScatterEvent & event) const override;

    virtual void LogInfo() const override;

    SpecularR(): BSDF(BXDFType(BSDF_REFLECTION | BSDF_SPECULAR)) { }

    virtual Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    // avoid calculate specular pdf directly
    Float Pdf(const SurfaceScatterEvent & event) const override;
};

class Dielectric : public  BSDF {
public:

    Dielectric(Float ior,bool enableT=true) : BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION |
                                                                            (enableT?BSDF_TRANSMISSION:0)) ),
                                                              ior(ior), enableT(enableT){
       invIor = 1.0f/ior;
    }

    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;

    Float eta( ) const override {
        return ior;
    }

private:
    Float  ior,invIor;
    bool  enableT;
};





class RoughDielectric : public  BSDF{
public:
    Spectrum f(const SurfaceScatterEvent & event) const override;

    Float Pdf(const SurfaceScatterEvent & event) const override;

    Spectrum sampleF(SurfaceScatterEvent & event, const vec2 & u) const override;

    void LogInfo( ) const override;
};

class Metal : public  BSDF{

};


class OrenNayar : public  BSDF{
public:
    OrenNayar(const Spectrum &R, Float sigma);

private:
    const Spectrum R;
    Float A, B;
};




