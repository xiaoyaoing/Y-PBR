#pragma  once

#include "Colors/Spectrum.hpp"
#include "Ray/Intersection.hpp"
#include "Ray/Ray.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"

#include <memory>
#include <optional>

class Scene;
class Primitive;
class Ray;

template<class T>
class Texture;


// LightFlags Declarations
enum LightFlags : int {
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


class Light {
public:
    Light(int _flags) : flags(_flags){};
    ///given a intersection,sample light
    virtual Spectrum sampleLi(const vec3 & ref, const vec2 &u,
                              vec3 *wi, Float *pdf,
                              Float * distance) const = 0;
    PositionAndDirectionSample sampleLi(const vec3 & ref, const vec2 &u) const{
        PositionAndDirectionSample result{};
        vec3 wi;Float pdf,distance;
        Spectrum  Li = sampleLi(ref,u,&wi,&pdf,&distance);
        if(isBlack(Li))
            return result;
        result.ray.o = ref + wi * distance;
        result.ray.d = -wi;
        result.posPdf = pdf;
        result.weight = Li/pdf;
        return result;
    }
    ///sample from light
    virtual  PositionAndDirectionSample sampleDirect(const vec2 & positionSample, const vec2 & dirSample) const = 0;
    virtual void pdfDirect(const Ray & ray,const vec3 & n,Float * pdfPos,Float * pdfDir) const  { };

    ///compute environment radiance
    virtual Spectrum Le(const Ray & ray) const  {return Spectrum(0);};
    /// returns total power of the light
    virtual Spectrum Power() { return Spectrum(); }
    virtual  void Preprocess(const Scene & scene) {

    }
    ///given a intersection compute lightPdf.Mainly used in mis
    virtual  Float PdfLi(const Intersection & pShape, const vec3 & ref) const {
        _NOT_IMPLEMENT_ERROR
    }
    virtual  bool isDeltaLight() const {return false;}

    virtual std::optional<Intersection> intersect(Ray & ray) const = 0;
    const int  flags;
};


class AreaLight : public  Light{
    Spectrum
    sampleLi(const vec3 & ref, const vec2 & u, vec3 * wi, Float * pdf, Float * distance) const override;

public:
    PositionAndDirectionSample sampleDirect(const vec2 & positionSample, const vec2 & dirSample) const override;

public:
    virtual  Spectrum directLighting(const Intersection & intr,const vec3 & wo) const ;


    AreaLight(const std::shared_ptr < Primitive > & _primitive,
                         const std::shared_ptr<Texture<Spectrum>> _emssision,
                         bool twoSide=false);
    Float PdfLi(const Intersection & pShape, const vec3 & ref) const override;

    void pdfDirect(const Ray &ray, const vec3 &n, Float *pdfPos, Float *pdfDir) const override;

    std::optional<Intersection> intersect(Ray & ray) const override;

protected:
    std::shared_ptr<Primitive> primitive;
    std::shared_ptr<Texture<Spectrum>>  emssision;
    bool twoSide;
    Float area;
};

class Infinite : public  Light{
public:
    std::optional < Intersection > intersect(Ray & ray) const override;
    void Preprocess(const Scene & scene) override;
    Infinite(): Light(int(LightFlags::Infinite)){}
    Infinite(LightFlags flag): Light(int(LightFlags::Infinite | flag)){}
protected:
    vec3 _worldCenter;
    Float _worldRadius;
};
