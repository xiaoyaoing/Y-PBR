#include "Light.hpp"
#include "../Primitives/Primitive.hpp"
#include "../scene.hpp"


#include "Common/Texture.hpp"
#include "Sampler/Warp.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"
#include "Primitives/Quad.hpp"

Spectrum
AreaLight::sampleLi(const vec3 & ref, const vec2 & u, vec3 * wi, Float * pdf, Float * distance) const {
    Intersection pShape = primitive->sample(ref, u, pdf, wi);
    if(*pdf==0){
        return Spectrum();
    }
    *distance = length(pShape.p - ref);
    return directLighting(pShape, -(*wi));

}



Spectrum AreaLight::directLighting(const Intersection & intr,const vec3 & wo) const {
    if( dot(intr.Ng, wo) <0){
        int k = 1;
    }
    return (twoSide || dot(intr.Ng, wo) > 0) ? emssision->eval(& intr) : Spectrum(0.f);
}



AreaLight::AreaLight(const std::shared_ptr < Primitive > & _primitive,
                     const std::shared_ptr<Texture<Spectrum>> _emssision,
                     bool twoSide)
                        :Light((int)LightFlags::Area),
                        primitive(_primitive),emssision(_emssision),twoSide(twoSide)
{

}

Float AreaLight::PdfLi(const Intersection & pShape, const vec3 & ref) const {
    return primitive->directPdf(pShape,ref);
}

PositionAndDirectionSample AreaLight::sampleDirect(const vec2 & positionSample, const vec2 & dirSample) const {
    PositionAndDirectionSample result;
    Intersection pshape = primitive->sample(positionSample, & result.posPdf);
    vec3 w;
    if(twoSide){
        vec2 u = dirSample;
        if(u[0]<0.5){
            u[0] *=2;
            w = Warp::squareToUniformHemisphere(u);
            result.dirPdf = Warp::squareToCosineHemispherePdf(w);
        }
        else {
            u[0] = u[0] *2-1;
            w = Warp::squareToUniformHemisphere(u);
            result.dirPdf = Warp::squareToCosineHemispherePdf(w);
            w=-w;
        }
    }
    else {
        w = Warp::squareToUniformHemisphere(dirSample);
        result.dirPdf = Warp::squareToCosineHemispherePdf(w);
    }
    result.n = pshape.Ng;

    vec3 s,t;
    coordinateSystem(result.n, s, t);
    result.ray = Ray(pshape.p,w.x * s+w.y * t+w.z*result.n);
    result.weight = directLighting(pshape, result.n);
    return result;
}

std::optional < Intersection > AreaLight::intersect(Ray & ray) const {
    return primitive->intersect(ray);
}

bool VisibilityTester::Unoccluded(const Scene & scene) const {
   // return true;
    vec3 dir =(p1.p-p0.p);
    Float distance = length(dir);
    dir= dir/distance;
    Ray ray(p0.p ,dir,p0.epsilon ,distance-p1.epsilon);
    return !scene.intersectP(ray);
}

std::optional < Intersection > Infinite::intersect(Ray & ray) const {
    Intersection its;
    its.p = ray.o + 2*_worldRadius;
    its.w = ray.d;
    return {its};
}

void Infinite::Preprocess(const Scene & scene) {
    scene.getWorldBound().BoundingSphere(& _worldCenter, & _worldRadius);
}
