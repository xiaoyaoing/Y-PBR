#include "Primitive.hpp"
class Sphere : public Primitive
{
public:
    Sphere(double radius, std::shared_ptr<Bsdf> bsdf);

    Intersection Sample(const vec2 & u, Float * pdf) const override;

    virtual std::optional<Intersection> intersect(Ray& ray) const ;

    virtual vec3 operator()(Float u, Float v) const;

    virtual vec3 normal(const vec3& pos) const;

    virtual void transform(const Transform &T);



protected:
    virtual void computeArea();
    virtual void computeBoundingBox();

private:
    vec3 origin;
    Float radius;
};