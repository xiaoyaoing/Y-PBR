//2022/7/17

#include "Primitive.hpp"

Quadric::Quadric(const nlohmann::json &j, std::shared_ptr<Bsdf> bsdf):Primitive(bsdf) {

}

bool  Quadric::intersect( Ray& ray, Intersection& intersection) const{
    return false;
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

