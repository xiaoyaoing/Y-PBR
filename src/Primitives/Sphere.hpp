#include "Primitive.hpp"
#include "Common/Transform.hpp"

class Sphere : public Primitive
{
public:
    Sphere(double radius, std::shared_ptr<BSDF> bsdf);
    Intersection sample(const vec2 & u, Float * pdf) const override;
    virtual std::optional<Intersection> intersect(Ray& ray) const ;
    virtual vec3 operator()(Float u, Float v) const;
    virtual vec3 normal(const vec3& pos) const;
    virtual void transform(const mat4 & T);

    bool occluded(const Ray & ray) const override;


protected:
    virtual void computeArea();
    virtual void computeBoundingBox();

    vec3 origin;
    Float radius;
};