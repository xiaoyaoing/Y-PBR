#pragma  once

#include "nlohmann/json.hpp"
#include "../Colors/Spectrum.hpp"
#include "BsdfLobes.hpp"


// BSDF Declarations
enum BXDFType {
    BSDF_REFLECTION = 1 << 0,
    BSDF_TRANSMISSION = 1 << 1,
    BSDF_DIFFUSE = 1 << 2,
    BSDF_GLOSSY = 1 << 3,
    BSDF_SPECULAR = 1 << 4,
    BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
               BSDF_TRANSMISSION,
};

class BXDF;

class Bsdf{

public:
    Spectrum f(const vec3 &wo , const vec3 &wi,
               BXDFType flags = BSDF_ALL) const ;

    ///sample a new dir
    ///wi wo in local space
    virtual Spectrum sampleF(const vec3 &wo, vec3 *wi, const vec2 &u,
                             Float *pdf, BXDFType type,BXDFType *sampledType) const ;

    int NumComponents(BXDFType flags = BSDF_ALL) const;

    void Add(BXDF *b) {
        BXDFs[nBXDFs++] = b;
    }



    void LogInfo();

    std::string name;

//    Float eta(){
//        return BXDF[]
//    }
private:
    BXDF *BXDFs[8];
    int nBXDFs = 0;
    //for debug
};

class BXDF{
public:

    BXDF(BXDFType type) :type(type) {};

    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const =0;

    virtual Float Pdf(const vec3 & wo,const vec3 & wi) const =0;

    virtual Spectrum sampleF(const vec3 &wo, vec3 *wi, const vec2 &u,
                             Float *pdf, BXDFType *sampledType) const =0;


    virtual void LogInfo() const =0;

    bool MatchesFlags(BXDFType t) const { return (type & t) == type; }

    BXDFType type;

};


class Mirror : public BXDF{

};

class LambertainR : public  BXDF{

public:

    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const override;

    Float Pdf(const vec3 & wo, const vec3 & wi) const override;

    virtual void LogInfo() const;

    LambertainR(Spectrum &  albedo);

    virtual Spectrum sampleF(const vec3 &wo, vec3 *wi, const vec2 &u,
                             Float *pdf, BXDFType *sampledType) const override;
private:
    bool useCosineSample;
    Spectrum albedo;
};

class LambertainT : public  BXDF{
public:

    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const override;

    virtual void LogInfo() const override;

    LambertainT(Spectrum &  albedo);

    virtual Spectrum sampleF(const vec3 &wo, vec3 *wi, const vec2 &u,
                             Float *pdf, BXDFType *sampledType) const override;

private:

    Spectrum albedo;
};

class SpecularR : public BXDF{
public:


    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const override;

    virtual void LogInfo() const override;

    SpecularR();

    virtual Spectrum sampleF(const vec3 &wo, vec3 *wi, const vec2 &u,
                             Float *pdf, BXDFType *sampledType) const override;

    // avoid calculate specular pdf directly
    Float Pdf(const vec3 & wo, const vec3 & wi) const override;
};

class Dielectric : public  BXDF {
public:

    Dielectric(Float ior,Spectrum albedo,bool enableT=true) : BXDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION |
                                                                (enableT?BSDF_TRANSMISSION:0)) ),
                                                        ior(ior),albedo(albedo),enableT(enableT){
       invIor = 1.0f/ior;
    }

    Spectrum f(const vec3 & wo, const vec3 & wi) const override;

    Float Pdf(const vec3 & wo, const vec3 & wi) const override;

    Spectrum sampleF(const vec3 & wo, vec3 * wi, const vec2 & u, Float * pdf, BXDFType * sampledType) const override;

    void LogInfo( ) const override;

private:
    Spectrum albedo;
    Float  ior,invIor;
    bool  enableT;
};

class Conductor : public BXDF {
public:
    Conductor();

    Spectrum f(const vec3 & wo, const vec3 & wi) const override;

    Float Pdf(const vec3 & wo, const vec3 & wi) const override;

    Spectrum sampleF(const vec3 & wo, vec3 * wi, const vec2 & u, Float * pdf, BXDFType * sampledType) const override;

    void LogInfo( ) const override;
};

class OrenNayar : public  BXDF{
public:
    OrenNayar(const Spectrum &R, Float sigma);

private:
    const Spectrum R;
    Float A, B;
};




