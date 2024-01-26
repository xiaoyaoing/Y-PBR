#include "Primitive.hpp"
#include "Common/Transform.hpp"

#pragma once

class Quad : public Primitive {
public:
    Quad(const Json& json) : Primitive(json) {}
    virtual std::optional<Intersection> intersect(Ray& ray) const;
    virtual vec3                        operator()(Float u, Float v) const {
        return _base + _edge0 * u + _edge1 * v;
    }
    virtual vec3 normal(const vec3& pos) const;
    void         transform(const mat4& T) override;
    Intersection sample(const vec2& u, Float* pdf) const override;

    Float directPdf(const Intersection& pShape, vec3 ref) const override;

    bool occluded(const Ray& ray) const override;

    Float powerToRadianceScale() const override;

protected:
    virtual void computeArea();
    virtual void computeBoundingBox();
    void         preCompute();

private:
    vec3 _base;
    vec3 _edge0, _edge1;
    vec2 _invSq;
    vec3 _normal;
};