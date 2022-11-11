#pragma  once

#include "../Colors/Spectrum.hpp"
#include "../Ray/Intersection.hpp"
#include "Ray/Ray.hpp"

class Scene;
class Primitive;
class Ray;

template<class T>
class Texture;

// LightFlags Declarations
enum class LightFlags : int {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8,
    Distant = 16
};

class VisibilityTester {
public:
    VisibilityTester() {}
    // VisibilityTester Public Methods
    VisibilityTester(const Intersection &p0, const Intersection &p1)
            : p0(p0), p1(p1) {}
    const Intersection &P0() const { return p0; }
    const Intersection &P1() const { return p1; }
    bool Unoccluded(const Scene  &scene) const;
//    Spectrum Tr(const Scene &scene, Sampler &sampler) const;

private:
   Intersection p0, p1;
};


struct LightSampleResult{
    vec3 lightN;
    Float lightPosPdf;
    Float lightDirPdf;
    Spectrum radiance;
    Ray ray;
};

class Light {
public:
    Light(int _flags) : flags(_flags){};
    ///given a intersection,sample light
    virtual Spectrum sampleLi(const Intersection &ref, const vec2 &u,
                              vec3 *wi, Float *pdf,
                              VisibilityTester *vis) const = 0;
    ///sample from light
    virtual  LightSampleResult sampleDirect(const vec2 & positionSample, const vec2 & dirSample) = 0;
    ///compute environment radiance
    virtual Spectrum Le(const Ray & ray) const  {return Spectrum(0);};
    /// returns total power of the light
    virtual Spectrum Power() { return Spectrum(); }
    virtual  void Preprocess(const Scene & scene) {}
    ///given a intersection compute lightPdf.Mainly used in mis
    virtual  Float PdfLi(const Intersection & pShape, const vec3 & ref) const {
        _NOT_IMPLEMENT_ERROR
    }
    virtual  bool isDeltaLight() const {return false;}
    const int  flags;
};


class AreaLight : public  Light{
    Spectrum
    sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;

public:
    LightSampleResult sampleDirect(const vec2 & positionSample, const vec2 & dirSample) override;

public:
    virtual  Spectrum directLighting(const Intersection & intr,const vec3 & wo) const ;


    AreaLight(const std::shared_ptr < Primitive > & _primitive,
                         const std::shared_ptr<Texture<Spectrum>> _emssision,
                         bool twoSide=false);
    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;


protected:
    std::shared_ptr<Primitive> primitive;
    std::shared_ptr<Texture<Spectrum>>  emssision;
    bool twoSide;
    Float area;
};
