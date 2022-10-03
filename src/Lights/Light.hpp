#pragma  once

#include "../Colors/Spectrum.hpp"
#include "../Ray/Intersection.hpp"

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
    Infinite = 8
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
    virtual Spectrum Sample_Li(const Intersection &ref, const vec2 &u,
                                vec3 *wi, Float *pdf,
                               VisibilityTester *vis) const = 0;

    Light(int _flags) : flags(_flags){};

    /// compute radiance in position with dir w
    /// \param intr info
    /// \return radiance
    virtual Spectrum directLighting(const Intersection & intr) const=0;

    ///compute environment radiance
    virtual Spectrum environmentLighting(const Ray & ray) const  {return Spectrum(0);};

    /// returns total power of the light
    virtual  Float Power() { return 0; }
    virtual  void Preprocess(const Scene & scene) {}
    virtual  Float directPdf(const Intersection & pShape, const vec3 & ref) const {
        _NOT_IMPLEMENT_ERROR
    }
//    virtual  Spectrum  directLighting(const Intersection & pShape) const {return Spectrum(0);}
    virtual  bool isDeltaLight() const {return false;}

    const int  flags;
};


class AreaLight : public  Light{
    Spectrum
    Sample_Li(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;

public:
    Spectrum directLighting(const Intersection & intr) const override;


    AreaLight(const std::shared_ptr < Primitive > & _primitive,
                         const std::shared_ptr<Texture<Spectrum>> _emssision,
                         bool twoSide=false);
    Float directPdf(const Intersection & pShape, const vec3 & ref) const override;


protected:
    std::shared_ptr<Primitive> primitive;
    std::shared_ptr<Texture<Spectrum>>  emssision;
    bool twoSide;
    Float area;
};
