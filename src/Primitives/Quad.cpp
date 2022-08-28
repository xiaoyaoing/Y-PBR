//2022/7/17

#include "Quad.hpp"

Quad::Quad(const nlohmann::json &j, std::shared_ptr<Bsdf> bsdf):Primitive(bsdf) {
}

void Quad::preCompute( ) {
    computeArea();
    computeBoundingBox();
    _invSq = vec2(1/ length2(_edge0),1/ length2(_edge1));
}

Quad::Quad(std::shared_ptr < Bsdf > bsdf):Primitive(bsdf) {
    _base = vec3 (0,0,0);
    _edge0 = vec3(1,0,0);
    _edge1 = vec3 (0,0,1);
    preCompute();
}

std::optional < Intersection > Quad::intersect(Ray & ray) const{
    if(bsdf->name=="ceiling"){
        int k=1;
    }
    Intersection its;
    vec3 n = normal(_base);
    Float dotW = dot(ray.d,n);
    if(std::abs(dotW)<Constant::EPSILON){
        return std::nullopt;
    }
    Float t = dot(n,(_base - ray.o))/dotW;

    if(bsdf->name=="light"){
        int k=1;
    }
    if (t < ray.nearT || t > ray.farT)
    return std::nullopt;

    vec3 q = ray(t);
    vec3 v = q - _base;
    Float l0 = dot(v,_edge0) * _invSq.x;
    Float l1 = dot(v,_edge1) * _invSq.y;


    if (l0 < 0.0f || l0 > 1.0f || l1 < 0.0f || l1 > 1.0f)
        return std::nullopt;

    ray.farT = t;
    its.p=q;
    its.setNormal(normal(its.p));
    its.primitive=this;
    its.bsdf=bsdf.get();
    return {its};
}

Intersection Quad::Sample(const vec2 & u, Float * pdf) const {
    Intersection its;
    its.p=_base + _edge0 * u[0] + _edge1 * u[1];
    its.setNormal(normal(its.p));
    *pdf = inv_area;
    return its;
}


vec3 Quad::normal(const vec3& pos) const{
    return normalize(cross(_edge1,_edge0));
}
void Quad::transform(const Transform &T) {
    if(bsdf->name=="leftWall"){
       int k=1;
    }
    std::string s = Mat4ToStr(T.matrix);
    _base = T * vec3(0.0f) ;
    _edge0 = T.transformVector(vec3(1.0f, 0.0f, 0.0f));
    _edge1 = T.transformVector(vec3(0.0f, 0.0f, 1.0f));
    _base -= _edge0*0.5f;
    _base -= _edge1*0.5f;

    preCompute();
}

void Quad::computeBoundingBox( ) {
    vec3 pMin = _base;
    vec3 pMax= _base+_edge1+_edge0;
    Bounds3 bounds;
    bounds = Union(bounds,pMin);
    bounds = Union(bounds,pMax);
    BB_ = bounds;
}

void Quad::computeArea() {
    area = glm::length(cross(_edge0 , _edge1));
    inv_area=1/area;
}


