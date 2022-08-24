//2022/7/17

#include "Sphere.hpp"
#include "../Sampler/Warp.hpp"

Sphere::Sphere(double radius, std::shared_ptr < Bsdf > bsdf) : Primitive(bsdf), radius(radius), origin(0.0) {
    computeArea();
    computeBoundingBox();
}

std::optional < Intersection > Sphere::intersect(Ray & ray) const {
    Intersection intersection;
    vec3 p = ray.o - origin;
    float B = dot(p, ray.d);
    float C = length2(p) - pow2(radius);
    float detSq = B * B - C;
    if ( detSq >= 0.0f ) {
        float det = std::sqrt(detSq);
        float t = - B - det;

        if ( t > ray.farT || t + 2 * det < ray.nearT ) {
            return std::nullopt;
        }

        if ( t < ray.nearT )
            t = t + 2 * det;

        ray.farT=t;
        intersection.p = ray(t);
        intersection.setNormal(normal(intersection.p));
        intersection.primitive = this;
        intersection.bsdf = bsdf.get();
        return {intersection};
    }
    return std::nullopt;

}

vec3 Sphere::operator ()(Float u, Float v) const {

}

vec3 Sphere::normal(const vec3 & pos) const {
    return normalize(pos - origin);
}

void Sphere::transform(const Transform & T) {
    origin = T * vec3(0,0,0);
    vec3  scale = extractScale(T.matrix) * vec4(vec3(1.0),0);
    radius = max(scale);
    computeBoundingBox();
    computeArea();

}

void Sphere::computeArea( ) {
    area = 4 * Constant::PI * radius * radius;
    inv_area = 1 / area;
}

void Sphere::computeBoundingBox( ) {
    BB_ = Bounds3(origin-radius,origin+radius);
    //todo
}

Intersection Sphere::Sample(const vec2 & u, Float * pdf) const {
    Intersection it;
    it.p = origin + radius * Warp::squareToUniformSphere(u);
    it.setNormal(normal(it.p));
    * pdf = inv_area;
    return it;
}


