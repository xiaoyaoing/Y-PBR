///**
// * @file
// * @author JunPing Yuan
// * @brief
// * @version 0.1
// * @date 2023/4/13
// *
// * @copyright Copyright (c) 2022
// *
// */
//
//#include "Reflection.hpp"
//#include "MicrofacetDistribution.hpp"
//
//template<class T>
//inline T rcp(T x){
//    return 1/x;
//}
//
//template<class T>
//inline T rsqrt(T x){
//    return 1/sqrt(x);
//}
//
//template<class T1,class T2>
//inline  auto select(bool mask,T1 a,T2 b){
//    return mask?a:b;
//}
//
//class Sampler;
//class RoughHair : public  BSDF{
//public:
//    RoughHair();
//
//    RoughHair(std::shared_ptr<Sampler> sampler);
//
//    Float Pdf(const SurfaceEvent & event) const override;
//
//protected:
//    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;
//    Spectrum f(const SurfaceEvent & event) const override;
//    Spectrum evalR(const vec3 & out,const vec3 & in) const;
//    Spectrum evalTTAndTRT(const vec3 & out,const vec3 & in) const;
//
//
//
//    Float sintheta(const vec3& w) const {
//        return w.y;
//    }
//
//    /* returns cos_theta */
//    Float costheta(const vec3& w) const {
//        return sqrt(sqr(w.x) + sqr(w.z));
//    }
//
//    /* returns tan_theta */
//    Float tantheta(const vec3& w) const {
//        return sintheta(w) / costheta(w);
//    }
//
//    /* extract theta coordinate from 3D direction
//     * -pi < theta < pi */
//    Float dir_theta(const vec3& w) const {
//        return atan2(sintheta(w), costheta(w));
//    }
//
//    /* extract phi coordinate from 3D direction.
//     * -pi < phi < pi
//     * Assuming phi(wi) = 0 */
//    Float dirPhi(const vec3& w) const {
//        return atan2(w.x, w.z);
//    }
//
//    Spectrum getSpectrum() const {
//        return Spectrum(612.0,549.0,465.0);
//    }
//
//    std::pair<Float,Float> sincos(Float  angle) const {
//        return {sin(angle), cos(angle)};
//    }
//
//    inline vec3 sph_dir(Float theta, Float gamma) const {
//        auto [sin_theta, cos_theta] = sincos(theta);
//        auto [sin_gamma,   cos_gamma]   = sincos(gamma);
//        return vec3 (sin_gamma * cos_theta, sin_theta, cos_gamma * cos_theta);
//    }
//
//    // NDF
//    Float D(const vec3 &m, const vec3 &h) const {
//        Float cos_theta = dot(h, m),
//                result;
//        result = m_roughness_squared * rcp(M_PI * sqr(1.f + (m_roughness_squared - 1.f) * sqr(cos_theta)));
//        return select(result * cos_theta > 1e-20f, result, 0.f);
//    }
//
//    Float smith_g1_visible(const vec3 &v, const vec3 &m, const vec3 &h) const {
//        Float cos_vm = dot(v, m),
//                result;
//        result = 2.f * rcp(cos_vm + sqrt(m_roughness_squared + (1.f - m_roughness_squared) * sqr(cos_vm)));
//
//        /* Assume consistent orientation (can't see the back
//           of the microfacet from the front and vice versa) */
//        if(dot(v, h) <= 0.f || cos_vm <= 0.f) result = 0.f;
//        return result;
//    }
//
//    Float smith_g1(const vec3 &v, const vec3 &m, const vec3 &h) const {
//        Float cos_vm = dot(v, m),
//                tmp, result;
//
//        result = 2.f * rcp(1.f + sqrt(m_roughness_squared * rcp(sqr(cos_vm)) + 1.f - m_roughness_squared));
//
//        /* Assume consistent orientation (can't see the back
//           of the microfacet from the front and vice versa) */
//        if(dot(v, h) <= 0.f || cos_vm <= 0.f)
//            result = 0.f;
//        return result;
//
//    }
//
//    Float G_(const vec3 &wi, const vec3 &wo, const vec3 &m, const vec3 &h) const {
//        return smith_g1_(wi, m, h) * smith_g1_(wo, m, h);
//    } Float G(const vec3 &wi, const vec3 &wo, const vec3 &m, const vec3 &h) const {
//        return smith_g1(wi, m, h) * smith_g1(wo, m, h);
//    }
//
//
//
//
//    Float smith_g1_(const vec3 &v, const vec3 &m, const vec3 &h) const {
//        return (dot(v, h) > 0 && dot(v, m) > 0);
//    }
//
//
//
//    MicrofacetDistribution * distr;
//    Float m_roughness,m_roughness_squared;
//    Float m_tan_tilt;
//    Float m_tilt;
//    Sampler * m_sampler;
//    vec2 alphaXY;
//    Float m_eta,m_inv_eta;
//    Spectrum sigmaA;
//    bool m_analytical;
//};