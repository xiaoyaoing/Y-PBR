#include "Disk.hpp"

std::optional < Intersection > Disk::intersect(Ray & ray) const {
    return Primitive::intersect(ray);
}

Float Disk::directPdf(const Intersection & pShape, vec3 ref) const {
    return Primitive::directPdf(pShape, ref);
}

void Disk::load(const Json & json, const Scene & scene) {
    Primitive::load(json, scene);
}

void Disk::computeArea( ) {
    Primitive::computeArea();
}

void Disk::computeBoundingBox( ) {
    Primitive::computeBoundingBox();
}

Intersection Disk::sample(const vec3 &ref, const vec2 &u, Float *pdf, vec3 *wi) const {
    return Primitive::sample(ref, u, pdf, wi);
}

Intersection Disk::sample(const vec2 & u, Float * pdf) const {
    return Primitive::sample(u, pdf);
}

Float Disk::powerToRadianceScale( ) const {
    return Primitive::powerToRadianceScale();
}
