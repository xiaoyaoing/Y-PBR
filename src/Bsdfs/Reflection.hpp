#pragma  once

#include "nlohmann/json.hpp"
#include "../Colors/Spectrum.hpp"
#include "BsdfLobes.hpp"


// BSDF Declarations
enum BxDFType {
    BSDF_REFLECTION = 1 << 0,
    BSDF_TRANSMISSION = 1 << 1,
    BSDF_DIFFUSE = 1 << 2,
    BSDF_GLOSSY = 1 << 3,
    BSDF_SPECULAR = 1 << 4,
    BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
               BSDF_TRANSMISSION,
};

class Bxdf;

class Bsdf{

public:
    Spectrum f(const vec3 &wo , const vec3 &wi,
               BxDFType flags = BSDF_ALL) const ;

    int NumComponents(BxDFType flags = BSDF_ALL) const;

    void Add(Bxdf *b) {
        bxdfs[nBxDFs++] = b;
    }
private:
    Bxdf *bxdfs[8];
    int nBxDFs = 0;
};

class Bxdf{
public:
    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const =0;

    virtual void LogInfo() const =0;

    bool MatchesFlags(BxDFType t) const { return (type & t) == type; }
    std::string name; //for debug
    BxDFType type;

};


class Mirror : public Bxdf{

};

class Lambertain : public  Bxdf{

public:

    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const override;

    virtual void LogInfo() const;

    Lambertain(Spectrum &  albedo);
private:

    Spectrum albedo;
};




