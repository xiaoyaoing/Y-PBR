#include "Primitive.hpp"
#pragma  once

class Quad : public Primitive
{
public:
    Quad(const nlohmann::json &j, std::shared_ptr<Bsdf> material);
    Quad(std::shared_ptr<Bsdf> material);
    virtual std::optional < Intersection > intersect(Ray & ray) const;
    virtual vec3 operator()(Float  u, Float  v) const {
        return _base + _edge0 * u + _edge1 * v;
    }
    virtual vec3 normal(const vec3& pos) const;
    virtual void transform(const Transform &T);
    Intersection Sample(const vec2 & u, Float * pdf) const override;

protected:
    virtual void computeArea();
    virtual void computeBoundingBox() ;
    void preCompute();

private:
    vec3  _base;
    vec3 _edge0, _edge1;
    vec2 _invSq;
};