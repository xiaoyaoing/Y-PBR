//2022/7/17

#include "Primitive.hpp"

Quadric::Quadric(const nlohmann::json &j, std::shared_ptr<Bsdf> bsdf):Primitive(bsdf) {

}

std::optional < Intersection > Quadric::intersect(Ray & ray) const{
    return std::nullopt;
}
vec3 Quadric::operator()(double u, double v) const{
    return vec3(0);
}
vec3 Quadric::normal(const vec3& pos) const{
    return  vec3(0);
}
void Quadric::transform(const Transform &T) {

}


void Quadric::computeArea() {

}


