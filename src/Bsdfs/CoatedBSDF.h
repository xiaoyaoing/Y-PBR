//
// Created by pc on 2023/10/30.
//

#pragma once
#include "Reflection.hpp"
#include "Fresnel.hpp"

//Coated 在任何表面材质加了一层dielectric 忽略光线在折射时发生的方向变化

template<class T>
class CoatedBSDF : public BSDF {
public:
    CoatedBSDF(std::shared_ptr<T> bsdf, std::shared_ptr<T> coating) : BSDF(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION)),
                                                                      bsdf(bsdf){};

    Float Pdf(const SurfaceEvent& event) const override;

    Float eta(const SurfaceEvent& event) const override;

protected:
    Spectrum sampleF(SurfaceEvent& event, const vec2& u) const override;

    Spectrum f(const SurfaceEvent& event) const override;

public:
    Float                 ior, invIor;
    std::shared_ptr<BSDF> bsdf;
};

template<class T>
Float CoatedBSDF<T>::Pdf(const SurfaceEvent& event) const {
    const vec3& wo = event.wo;
    const vec3& wi = event.wi;

    auto F = Fresnel::dielectricReflectance(wo.z < 0 ? ior : 1.f / ior, abs(wo.z));

    return bsdf->Pdf(event) * (SameHemisphere(wo, wi) ? F : (1 - F));
}

template<class T>
Float CoatedBSDF<T>::eta(const SurfaceEvent& event) const {
    if (event.wi.z * event.wo.z > 0)
        return 1;
    return event.wo.z < 0 ? ior : invIor;
}

template<class T>
Spectrum CoatedBSDF<T>::sampleF(SurfaceEvent& event, const vec2& u) const {
    const vec3& wo = event.wo;
    auto        F  = Fresnel::dielectricReflectance(wo.z < 0 ? ior : 1.f / ior, abs(wo.z));
    if (u.x < F) {
        event.wi         = Reflect(wo, vec3(0, 0, 1));
        event.sampleType = BXDFType(BSDF_SPECULAR | BSDF_REFLECTION);
        event.pdf        = F;
        return Spectrum(F);
    } else {
        vec2     remapU = vec2((u.x - F) / (1 - F), u.y);
        Spectrum s      = (1 - F) * bsdf->sampleF(event, remapU, true);
        event.pdf *= (1 - F);
        return s;
    }
    return Spectrum();
}

template<class T>
Spectrum CoatedBSDF<T>::f(const SurfaceEvent& event) const {
    const vec3& wo = event.wo;
    const vec3& wi = event.wi;

    auto F = Fresnel::dielectricReflectance(wo.z < 0 ? ior : 1.f / ior, abs(wo.z));

    return bsdf->f(event, true) * (SameHemisphere(wo,wi)?F:(1-F));

}