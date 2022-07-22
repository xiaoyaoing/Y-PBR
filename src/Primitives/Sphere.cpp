//2022/7/17

#include "Primitive.hpp"

Sphere::Sphere(double radius, std::shared_ptr<Bsdf> bsdf):Primitive(bsdf),radius(radius),origin(0.0){

}

bool  Sphere::intersect(Ray& ray, Intersection& intersection) const{
    vec3  p = ray.start - origin;
    glm::vec3 v=p;
    float B = dot(p,ray.direction);
    float C = length2(p) - pow2(radius);
    float detSq = B*B - C;
    if (detSq >= 0.0f) {
        float det = std::sqrt(detSq);
        float t = -B - det;
        if (t < ray.farT && t > ray.nearT) {
            ray.farT=t;
            intersection.pos=ray(t);
            intersection.normal= normal(intersection.pos);
//            data.primitive = this;
//            data.as<SphereIntersection>()->backSide = false;
            return true;
        }
        t = -B + det;
        if (t < ray.farT && t > ray.nearT) {
            ray.farT = t;
            intersection.pos=ray(t);
            intersection.normal= normal(intersection.pos);
//            data.primitive = this;
//            data.as<SphereIntersection>()->backSide = true;
            return true;
        }
    }

    return false;



    return false;
}
vec3 Sphere::operator()(double u, double v) const{

}
vec3 Sphere::normal(const vec3& pos) const{
    return  normalize(pos-origin);
}
void Sphere::transform(const Transform &T) {

}

void Sphere::computeArea() {

}

void Sphere::computeBoundingBox() {

}




