#pragma once
#include "Ray/Intersection.hpp"
#include "scene.hpp"
#include "Reflection.hpp"

Float FresnelMoment1(Float invEta);
Float FresnelMoment2(Float invEta);

class BSSRDF{
public:
    virtual Spectrum sampleS(const Scene & scene,Float u1,vec2 u2,const SurfaceEvent & po,Intersection * pi,
                             Float * pdf) const = 0;
    BSSRDF(Float eta) : eta(eta){}

    ///Used to determine whether the sampling point is on an object
   // void setBSDF(const BSDF * _bsdf) {bsdf = _bsdf;}
protected:
//    const BSDF * bsdf;
    Float eta;
    const Intersection * its;
};

///SeparableBSSRDF
class SeparableBSSRDF : public BSSRDF {
    friend  class  BSSRDFAdapater;
public:
    SeparableBSSRDF(Float eta) : BSSRDF(eta){}
    Spectrum sampleS(const Scene & scene,Float u1,vec2 u2,const SurfaceEvent & po,Intersection * pi,
                     Float * pdf) const override;
    Spectrum sampleSp(const Scene &scene, Float u1, const vec2 &u2,const SurfaceEvent & po,Intersection * pi,
                       Float *pdf) const;
    Float pdfSP(const SurfaceEvent &po, const Intersection & pi) const;

    Spectrum Sw(const vec3 &w) const ;

    // SeparableBSSRDF Interface
    virtual Spectrum Sr(Float d, const SurfaceEvent & event) const = 0;
    virtual Float sampleSr(int ch, Float u, const SurfaceEvent & event) const = 0;
    virtual Float pdfSr(int ch, Float r, const SurfaceEvent & event) const = 0;
    virtual Float       maxSr( int ch ,const SurfaceEvent & event) const = 0;
protected:
};

class DisneyBSSRDF : public  SeparableBSSRDF {
    Spectrum Sr(Float d, const SurfaceEvent & event) const override;

    Float sampleSr(int ch, Float u, const SurfaceEvent & event) const override;

    Float pdfSr(int ch, Float r, const SurfaceEvent & event) const override;

public:
    Float maxSr(int ch, const SurfaceEvent & event) const override;

private:

    std::shared_ptr<Texture<Spectrum>> scatterDistance,color;
public:
    DisneyBSSRDF(std::shared_ptr<Texture<Spectrum>> scatterDistance,
                 std::shared_ptr<Texture<Spectrum>> color,
                 Float eta):
                 SeparableBSSRDF(eta) ,scatterDistance(scatterDistance),color(color){}
};

struct BSSRDFTable{
    const int nRhoSamples, nRadiusSamples;
    std::unique_ptr<Float[]> rhoSamples, radiusSamples;
    std::unique_ptr<Float[]> profile;
    std::unique_ptr<Float[]> rhoEff;
    std::unique_ptr<Float[]> profileCDF;

    // BSSRDFTable Public Methods
    BSSRDFTable(int nRhoSamples, int nRadiusSamples);
    inline Float EvalProfile(int rhoIndex, int radiusIndex) const {
        return profile[rhoIndex * nRadiusSamples + radiusIndex];
    }
};

class TabelBSSRDF : public  SeparableBSSRDF {
public:
    Spectrum Sr(Float d, const SurfaceEvent & event) const override;

    Float sampleSr(int ch, Float u, const SurfaceEvent & event) const override;

    Float pdfSr(int ch, Float r, const SurfaceEvent & event) const override;

    Float maxSr(int ch, const SurfaceEvent & event) const override;

    TabelBSSRDF(Float eta, std::shared_ptr<Texture<Spectrum>> sigma_a,
            std::shared_ptr<Texture<Spectrum>> sigma_s, Float g);

private:
    std::shared_ptr<Texture<Spectrum>> sigmaA,sigmaS;
    BSSRDFTable  table;
    Float g;
};


class BSSRDFAdapater : public  BSDF{
public:
    BSSRDFAdapater(const SeparableBSSRDF * bssrdf) : bssrdf(bssrdf), BSDF(BXDFType(BSDF_DIFFUSE | BSDF_REFLECTION)){}

    Float Pdf(const SurfaceEvent & event) const override;

    Float eta(const SurfaceEvent & event) const override;

protected:
    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    Spectrum f(const SurfaceEvent & event) const override;

    const SeparableBSSRDF * bssrdf;
};


static Float BeamDiffusionMultScattering(Float sigmaS,Float sigmaA,Float g,Float eta,Float r);
static Float BeamDiffusionSingleScattering(Float sigmaS,Float sigmaA,Float g,Float eta,Float r);
static void computeBeamDiffusionBSSRDF(Float g,Float eta,BSSRDFTable * table);