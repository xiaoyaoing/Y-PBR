//2022/7/17

#include "Sphere.hpp"

Sphere::Sphere(double radius, std::shared_ptr<Bsdf> bsdf):Primitive(bsdf),radius(radius),origin(0.0){
    computeArea();
    computeBoundingBox();
}

std::optional<Intersection>  Sphere::intersect(Ray& ray) const{
    Intersection intersection;
    vec3  p = ray.o - origin;
    float B = dot(p,ray.d);
    float C = length2(p) - pow2(radius);
    float detSq = B*B - C;
    if (detSq >= 0.0f) {
        float det = std::sqrt(detSq);
        float t = -B - det;

        if(t>ray.farT || t+2*det<ray.nearT){
            return std::nullopt;
        }

        if(t<ray.nearT)
            t=t+2*det;

        intersection.p=ray(t);
        intersection.setNormal(normal(intersection.p));
        intersection.primitive=this;
        intersection.bsdf=bsdf.get();
            return {intersection};
        }


    return std::nullopt;

}
vec3 Sphere::operator()(double u, double v) const{

}
vec3 Sphere::normal(const vec3& pos) const{
    return  normalize(pos-origin);
}
void Sphere::transform(const Transform &T) {
    origin = T.position;
    radius = radius * ((T.scale.x + T.scale.y + T.scale.z) / 3.0);

}

void Sphere::computeArea() {
    area = Constant::PI * radius * radius;
    inv_area =  1/area;
}

void Sphere::computeBoundingBox() {
    //todo
}

Intersection Sphere::Sample(const vec2 & u, Float * pdf) const {
    Intersection it;

    return it;
}
