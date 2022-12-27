#include "Primitive.hpp"
#include "Common/Transform.hpp"

class Disk : public  Primitive {
public:
    std::optional < Intersection > intersect(Ray & ray) const override;

    Float directPdf(const Intersection & pShape, vec3 ref) const override;

    void load(const Json & json, const Scene & scene) override;

    Intersection sample(const vec3 & ref, const vec2 & u, Float * pdf) const override;

    Intersection sample(const vec2 & u, Float * pdf) const override;

    Float powerToRadianceScale( ) const override;

protected:
    void computeArea( ) override;

    void computeBoundingBox( ) override;
};