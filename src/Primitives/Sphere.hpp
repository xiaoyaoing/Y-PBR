#include "Primitive.hpp"
#include "Common/Transform.hpp"

class Sphere : public Primitive
{
public:
    Sphere(double radius, std::shared_ptr<BSDF> bsdf);
    Intersection sample(const vec2 & u, Float * pdf) const override;

    Intersection sample(const Intersection & ref, const vec2 & u, Float * pdf) const override;

    std::optional<Intersection> intersect(Ray& ray) const override;
    vec3 operator()(Float u, Float v) const override;
    vec3 normal(const vec3& pos) const override;
    void transform(const mat4 & T) override;
    bool occluded(const Ray & ray) const override;

    Float directPdf(const Intersection & pShape, vec3 ref) const override;

protected:
    virtual void computeArea();
    virtual void computeBoundingBox();

    vec3 center;
    Float radius;
};