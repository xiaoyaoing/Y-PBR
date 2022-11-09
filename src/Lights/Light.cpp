#include "Light.hpp"
#include "../Primitives/Primitive.hpp"
#include "../scene.hpp"


#include "Common/Texture.hpp"
#include "Sampler/Warp.hpp"
Spectrum
AreaLight::sampleLi(const Intersection & ref, const vec2 & u, vec3 * wi, Float * pdf, VisibilityTester * vis) const {

    Intersection pShape = primitive->sample(ref, u, pdf);

    if(*pdf==0){
        return Spectrum();
    }

    *wi = normalize(pShape.p-ref.p);
    if( hasNan(*wi))
        return Spectrum();
    auto t = pShape.p - ref.p;
    t = normalize(t);
    pShape.w = t;
    *vis = VisibilityTester(ref, pShape);
    return directLighting(pShape);

}



Spectrum AreaLight::directLighting(const Intersection & intr) const {
    return (twoSide || dot(intr.Ng, -intr.w) > 0) ? emssision->Evaluate(&intr): Spectrum(0.f);
}



AreaLight::AreaLight(const std::shared_ptr < Primitive > & _primitive,
                     const std::shared_ptr<Texture<Spectrum>> _emssision,
                     bool twoSide)
                        :Light((int)LightFlags::Area),
                        primitive(_primitive),emssision(_emssision),twoSide(twoSide)
{

}

Float AreaLight::directPdf(const Intersection & pShape, const vec3 & ref) const {
    return primitive->directPdf(pShape,ref);
}

LightSampleResult AreaLight::sampleDirect(const vec2 & positionSample, const vec2 & dirSample) {
    LightSampleResult result;
    Intersection pshape = primitive->sample(positionSample, & result.lightPosPdf);
    vec3 w;
    if(twoSide){
        vec2 u = dirSample;
        if(u[0]<0.5){
            u[0] *=2;
            w = Warp::squareToUniformHemisphere(u);
            result.lightDirPdf = Warp::squareToCosineHemispherePdf(w);
        }
        else {
            u[0] = u[0] *2-1;
            w = Warp::squareToUniformHemisphere(u);
            result.lightDirPdf = Warp::squareToCosineHemispherePdf(w);
            w=-w;
        }
    }
    else {
        w = Warp::squareToUniformHemisphere(dirSample);
        result.lightDirPdf = Warp::squareToCosineHemispherePdf(w);
    }
    result.lightN = pshape.Ng;
    result.radiance = directLighting(pshape);

    vec3 s,t;
    coordinateSystem(result.lightN,s,t);
    result.ray = Ray(pshape.p,w.x * s+w.y * t+w.z*result.lightN);
    return result;
}

bool VisibilityTester::Unoccluded(const Scene & scene) const {
    vec3 dir =(p1.p-p0.p);
    Float distance = length(dir);
    dir= normalize(dir);
    Ray ray(p0.p ,dir,Constant::EPSILON,distance-Constant::EPSILON);
    auto t = scene.intersect(ray);
   // return !t;
    return !scene.intersectP(ray);
}
