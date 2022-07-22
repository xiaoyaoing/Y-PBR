//2022/7/17

#include "Primitive.hpp"

Triangle::Triangle(const vec3& v0, const vec3& v1, const vec3& v2, std::shared_ptr<Bsdf> Bsdf) : Primitive(Bsdf){

}

Triangle::Triangle(const vec3& v0, const vec3& v1, const vec3& v2,
         const vec3& n0, const vec3& n1, const vec3& n2, std::shared_ptr<Bsdf> Bsdf) : Primitive(Bsdf){

}

bool  Triangle::intersect( Ray& ray, Intersection& intersection) const{
    return false;
}
vec3 Triangle::operator()(double u, double v) const{
    return vec3(0);
}
vec3 Triangle::normal(const vec3& pos) const{
    return vec3(0);
}
void Triangle::transform(const Transform &T) {

}


vec3 Triangle::interpolatedNormal(const glm::dvec2 &uv) const {
    return Primitive::interpolatedNormal(uv);
}

void Triangle::computeArea() {

}

void Triangle::computeBoundingBox() {

}

