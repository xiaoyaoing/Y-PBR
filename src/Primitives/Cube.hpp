#include "Primitive.hpp"

class Cube : public  Primitive{
public:
    Cube(std::shared_ptr<BSDF> material);
    virtual std::optional < Intersection > intersect(Ray & ray) const;
    virtual vec3 operator()(Float u, Float v) const;
    virtual vec3 normal(const vec3& pos) const;
    virtual void transform(const Transform &T);
    Intersection Sample(const vec2 & u, Float * pdf) const override;

    bool occluded(const Ray & ray) const override;

protected:
    virtual void computeArea();
    virtual void computeBoundingBox() ;
    void preCompute(){
        computeArea();
        computeBoundingBox();
    }

private:
    mat4 _rot;
    mat4 _invRot;
    vec3 _pos;
    vec3 _scale;
};




