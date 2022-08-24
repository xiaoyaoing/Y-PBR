#pragma  once

#include "../Colors/Spectrum.hpp"
#include "../Ray/Intersection.hpp"

class Scene;
class Primitive;

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


    /// compute radiance in position with dir w
    /// \param intr info
    /// \param w  dir
    /// \return radiance
    virtual Spectrum L(const Intersection &intr, const vec3 &w) const=0;
};


class AreaLight : public  Light{
    Spectrum
    Sample_Li(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const override;

public:
    Spectrum L(const Intersection & intr, const vec3 & w) const override;


    AreaLight(const std::shared_ptr < Primitive > & primitive,
                         const Spectrum & albedo,
                         bool twoSide=false);

private:
    std::shared_ptr<Primitive> primitive;
    Spectrum  albedo;
    bool twoSide;
    Float area;
};
