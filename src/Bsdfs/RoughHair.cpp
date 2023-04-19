/**
 * @file 
 * @author JunPing Yuan
 * @brief 
 * @version 0.1
 * @date 2023/4/13
 *
 * @copyright Copyright (c) 2022
 *
 */
#include "RoughHair.h"
#include "Common/Frame.hpp"
#include "Fresnel.hpp"
#include "Sampler/Sampler.hpp"
/* sample microfacets from a tilted mesonormal */
std::pair<vec3, Float> sample_wh(const vec3 &wi, const vec3 &wm,
                                     const MicrofacetDistribution * distr,
                                     const vec2 &sample1,const vec2 & alphaXY) {
    /* Coordinate transformation for microfacet sampling */
    Frame wmFrame;
    wmFrame.n = wm;
    wmFrame.bitTangent = cross(vec3(0.f, 1.f, 0.f), wm);
    wmFrame.tangent = cross(wmFrame.n, wmFrame.bitTangent);
    vec3 wh = distr->Sample_wh(wmFrame.toLocal(wi), sample1, alphaXY);
    auto pdf = distr->Pdf(wi, wh, alphaXY);
    return {wmFrame.toWorld(wh), pdf};

}


std::tuple<Float, Float, Float, Float> fresnel(Float cos_theta_i, Float eta) {
    bool outside_mask = cos_theta_i > 0;
    Float  rcp_eta = rcp(eta),
            eta_it = select(outside_mask, eta, rcp_eta),
            eta_ti = select(outside_mask, rcp_eta, eta);
    Float  cosThetaT;
    auto r = Fresnel::dielectricReflectance(1/eta,cos_theta_i,cosThetaT);
    return {r,cosThetaT,eta_it,eta_ti};
}

static inline vec3 refract(const vec3 &out, const vec3 &wh, Float cosThetaT,
                         Float eta) {
    auto whDotOut = dot(out,wh);
    return  (eta * whDotOut - (whDotOut>0?1:-1)*cosThetaT)* wh - eta* out ;
}

static inline  vec3  reflect(const vec3 & out,const vec3 & wh){
    return -out + 2 * dot(wh,out) * wh;
}



Float RoughHair::Pdf(const SurfaceEvent & event) const {
    return 1;
    vec3  wo = event.wi;
    vec3  wi = event.wo;
    // CAUTIOUS: this is only an estimation of the PDF, is now disabled
    //return 1.f;
//    std::runtime_error("error");
//    assert(false);
    // check visibility because of scale tilt
    Float phi_o = dirPhi(wo);
    /* dot(wi, wmi) > 0 */
    Float phi_m_max = acos(std::max(-m_tan_tilt * tantheta(wi), 0.f));
    if (isnan(phi_m_max))
        return 0.f;
    Float phi_m_min = -phi_m_max;

    /* dot(wo, wm) > 0 */
    Float tmp1 = acos(std::max(-m_tan_tilt * tantheta(wo), 0.f));
    if (isnan(tmp1))
        return 0.f;

    /* absorption */
//    Spectrum wavelengths = getSpectrum();
//    Spectrum mu_a = fmadd(m_pheomelanin, pheomelanin(wavelengths),
//                          m_eumelanin * eumelanin(wavelengths));

    /* dot(wo, wmi) > 0 */
    Float phi_m_max_r = tmp1 + phi_o;
    Float phi_m_min_r = -tmp1 + phi_o;

    vec3 wh = normalize(wi + wo);

    Float pdf_r(0.f), pdf_tt(0.f), pdf_trt(0.f);

    /* Construct a microfacet distribution matching the
       roughness values at the current surface position. */
//    MicrofacetDistribution distr(m_type, m_roughness, false);

    /* initial sample resolution */
    Float res = m_roughness * .8f;
    Float scale = (phi_m_max - phi_m_min) * .5f;
    size_t intervals = 2 * ceil(scale / res) + 1;
    /* modified resolution based on integral domain */
    res = (phi_m_max - phi_m_min) / Float(intervals);
    // integrate using Simpson's rule
    for (size_t i = 0; i < intervals; i++) {
        Float phi_mi = phi_m_min + i * res;
        vec3 wmi = sph_dir(m_tilt, phi_mi);
        /* R */
        /* sample wh1 */
        vec2 sample1 = m_sampler->getNext2D();
        auto [wh1, pdfh1] = sample_wh(wi, wmi, distr, sample1,alphaXY);
        auto [R1, cos_theta_t1, eta_it1, eta_ti1] = fresnel(dot(wi, wh1), Float(m_eta));
        /* TT */
        Float T1 = 1.f - R1;
        vec3 wt = refract(wi, wh1, cos_theta_t1, eta_ti1);
        Float phi_t = dirPhi(wt);
        Float phi_mt = 2.f * phi_t - phi_mi;
        vec3 wmt = sph_dir(-m_tilt, phi_mt);
        /* sample wh2 */
        vec2 sample2 = m_sampler->getNext2D();
        auto [wh2, pdfh2] = sample_wh(-wt, wmt, distr, sample2,alphaXY);
        Float R2 = std::get<0>(fresnel(dot(wh2, -wt), Float(m_inv_eta)));
        Spectrum At = exp(sigmaA * 2.f * cos(phi_t - phi_mi) / costheta(wt));
        Spectrum TT = T1 * (1.f - R2) * At;
        /* TRT */
        vec3 wtr = reflect(wt, wh2);
        Float phi_tr = dirPhi(wtr);
        Float twottrpi = -2.f * (phi_t - phi_tr) + M_PI;
        Float phi_mtr = phi_mi + twottrpi;
        vec3 wmtr = sph_dir(-m_tilt, phi_mtr);
        /* sample wh3 */
        vec2 sample3 = m_sampler->getNext2D();
        auto [wh3, pdfh3] = sample_wh(wtr, wmtr, distr, sample3,alphaXY);
        Float T3 = 1.f - std::get<0>(fresnel(dot(wh3, wtr), Float(m_inv_eta)));
        Spectrum Atr = exp(sigmaA * 2.f * cos(phi_tr - phi_mt) / costheta(wtr));
        Spectrum TRT = T1 * R2 * T3 * Atr;

        Float total_energy = luminace(Spectrum(R1)) + luminace(TT) + luminace(TRT);

        if (total_energy <= 0)
            continue;

        Float weight = (i == 0 || i == intervals - 1) ? 0.5f : (i % 2 + 1); /* Simpson */
        weight *= 0.5f * cos(phi_mi); /* cos(phi)dphi = dh, diameter = 2 */

        // R
        pdf_r += select(phi_mi < phi_m_max_r && phi_mi > phi_m_min_r,
                        std::max(0.f, R1 / total_energy * weight
                               * D(wmi, wh) * smith_g1_visible(wi, wmi, wh) * 0.25f),
                        0.f);
        // TT
        if (dot(wo, wt) > m_inv_eta) {
            vec3 wh2 = -wt + m_inv_eta * wo;
            Float rcp_norm_wh2 = rcp(wh2.length());
            wh2 *= rcp_norm_wh2;
            Float pdf_h2 = D(wmt, wh2) * smith_g1_visible(-wt, wmt, wh2) * dot(-wt, wh2);
            Float dwh2_wtt = sqr(m_inv_eta * rcp_norm_wh2) * dot(-wo, wh2);
            vec3 wmt_ = sph_dir(0, phi_mt);
            Float result = luminace(TT) / total_energy * pdf_h2 * dwh2_wtt * weight * G_(-wt, -wo, wmt_, wh2);
            if(isinf(result)) result = 0;
           // masked(result, !enoki::isfinite(result)) = 0;
            pdf_tt += result;
        }
        // TRT
        if (dot(-wtr, wo) > m_inv_eta) {
            vec3 wmtr_ = sph_dir(0, phi_mtr);
            vec3 wh3 = wtr + m_inv_eta * wo;
            Float rcp_norm_wh3 = rcp(wh3.length());
            wh3 *= rcp_norm_wh3;
            Float pdf_h3 = D(wmtr, wh3) * smith_g1_visible(wtr, wmtr, wh3) * dot(wtr, wh3);
            Float dwh3_wtrt = sqr(m_inv_eta * rcp_norm_wh3) * dot(-wo, wh3);
            Float result = luminace(TRT) / total_energy * pdf_h3 * dwh3_wtrt * weight * G_(wtr, -wo, wmtr_, wh3);
            if(isinf(result)) result = 0;

            //  masked(result, !enoki::isfinite(result)) = 0.f;
            pdf_trt += result;
        }

        return  (pdf_r + pdf_tt + pdf_trt) * (2.f * res / 3.f);

    }}

    Spectrum RoughHair::evalR(const vec3 &out, const vec3 &in) const {
        vec3 wi = out;
        vec3 wo = in;
        Spectrum R = Spectrum();

        vec3 wh = normalize(wi + wo);

        Float phi_o = dirPhi(wo);
        Float phi_h = dirPhi(wh);

        // compute valid phi_mi
        /* dot(wi, wmi) > 0 */
        Float phi_m_max1 = acos(std::max(-m_tan_tilt * tantheta(wi), 0.f));

        if (isnan(phi_m_max1))
            return Spectrum();
        Float phi_m_min1 = -phi_m_max1;

        /* dot(wo, wmi) > 0 */
        Float phi_m_max2 = acos(std::max(-m_tan_tilt * tantheta(wo), 0.f)) + phi_o;
        if (isnan(phi_m_max2))
            return Spectrum();
        Float phi_m_min2 = -phi_m_max2 + 2.f * phi_o;

        Float phi_m_min = std::max(phi_m_min1, phi_m_min2) + 1e-5f;
        Float phi_m_max = std::min(phi_m_max1, phi_m_max2) - 1e-5f;

        if (phi_m_min > phi_m_max)
            return Spectrum();

        Float integral = 0.f;
        Float d_max = phi_h - phi_m_max;
        Float d_min = phi_h - phi_m_min;
        if (m_analytical) { // TODO: beckmann
            if (m_tilt == 0.f) {
                Float A = (m_roughness_squared - 1) * sqr(costheta(wh));
                Float temp1 = A * rcp(A + 1) * (sin(2 * d_max) * rcp(A * cos(2 * d_max) + A + 2.f) -
                                                 sin(2 * d_min) * rcp(A * cos(2 * d_min) + A + 2.f));
                Float temp2 = (A + 2) * pow(A + 1, Float(-1.5)) *
                               (atan(tan(d_min) * rsqrt(A + 1)) - atan(tan(d_max) * rsqrt(A + 1)));
                integral = temp1 + temp2;
            } else {
                auto [sm, cm] = sincos(m_tilt);
                Float C = sqrt(1.f - m_roughness_squared);
                Float A = cm * costheta(wh) * C;
                Float B = sm * sintheta(wh) * C;
                Float A2 = sqr(A);
                Float B2 = sqr(B);
                Float tmp1 = rsqrt(sqr(B - 1.f) - A2);
                Float tmp2 = rsqrt(sqr(B + 1.f) - A2);

                auto [smax, cmax] = sincos(d_max);
                auto [smin, cmin] = sincos(d_min);
                Float tmax = smax / (1.f + cmax);
                Float tmin = smin / (1.f + cmin);

                Float temp1 = 2.f * (A2 - B2 + 3.f * B - 2) * sqr(tmp1) * tmp1 *
                               (atan((A - B + 1.f) * tmp1 * tmax) -
                                atan((A - B + 1.f) * tmp1 * tmin));
                Float temp2 = 2.f * (A2 - B2 - 3.f * B - 2) * sqr(tmp2) * tmp2 *
                               (atan((B - A + 1.f) * tmp2 * tmax) -
                                atan((B - A + 1.f) * tmp2 * tmin));
                Float temp3 = A * sqr(tmp1) *
                               (smax / (A * cmax + B - 1.f) -
                                smin / (A * cmin + B - 1.f));
                Float temp4 = A * sqr(tmp2) *
                               (smax / (A * cmax + B + 1.f) -
                                smin / (A * cmin + B + 1.f));

                integral = 0.5f * (temp1 + temp2 + temp3 + temp4);
            }
            integral *= m_roughness_squared * Constant::INV_TWO_PI;
        } else { /* falls back to numerical integration */
            /* initial sample resolution */
            Float res = m_roughness * .7f;
            Float scale = (phi_m_max - phi_m_min) * .5f;
            size_t intervals = 2 * ceil(scale / res);
            /* modified resolution based on integral domain */
            res = (phi_m_max - phi_m_min) / Float(intervals);
            // integrate using Simpson's rule
            for (size_t i = 0; i <= intervals; i++) {
                Float phi_m = phi_m_min + i * res;
                vec3 wm = sph_dir(m_tilt, phi_m);
                Float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);
                integral += weight * D(wm, wh) * G(wi, wo, wm, wh) * G_(wi, wo, vec3(wm.x, 0.f, wm.z), wh);
            }
            integral *= (2.f / 3.f * res);
        }
        auto temp = dot(wi, wh);
        Float F = std::get<0>(fresnel(dot(wi, wh), Float(m_eta)));
        R = Spectrum (0.25f * F * std::max(0.f, integral));

        return R; 
    }
    Spectrum RoughHair::evalTTAndTRT(const vec3 &out, const vec3 &in) const {

        auto wi = out;
        auto wo = in;
        /* dot(wi, wmi) > 0 */
        Float phi_m_max = acos(std::max(-m_tan_tilt * tantheta(wi), 0.f));
        if (isnan(phi_m_max))
            return Spectrum();
        Float phi_m_min = -phi_m_max;

        /* dot(wo, wmo) < 0 */
        Float tmp1 = acos(std::min(m_tan_tilt * tantheta(wo), 0.f)); //x
        if (isnan(tmp1))
            return Spectrum();

        Float res = m_roughness * .8f;

//        /* absorption */
//        Spectrum wavelengths = get_spectrum(si);
//        Spectrum mu_a = fmadd(m_pheomelanin, pheomelanin(wavelengths),
//                              m_eumelanin * eumelanin(wavelengths));

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position.
        */
      //  MicrofacetDistribution distr(m_type, m_roughness, true);
//        if (m_type == MicrofacetType::Beckmann) {
//            /* sample_visible = true would be too slow for beckmann */
//            distr = MicrofacetDistribution(m_type, m_roughness, false);
//        }

        Float scale = (phi_m_max - phi_m_min) * .5f;
        size_t intervals = 2 * ceil(scale / res);
        res = (phi_m_max - phi_m_min) / intervals;
        Spectrum S_tt = Spectrum(0), S_trt =Spectrum(0);
        for (size_t i = 0; i <= intervals; i++) {
            Float phi_mi = phi_m_min + i * res;
            vec3 wmi = sph_dir(m_tilt, phi_mi);

            /* sample wh1 */
            vec2 sample1 = m_sampler->getNext2D();
            vec3 wh1 = std::get<0>(sample_wh(wi, wmi, distr, sample1,alphaXY));
          //  wh1 = vec3 (-0.65917182, -0.f858450755, 0.747076332);
            Float cos_ih1 = dot(wi, wh1);
            if (!(cos_ih1 > 1e-5f))
                continue;

            /* fresnel coefficient */
            auto [R1, cos_theta_t1, eta_it1, eta_ti1] = fresnel(dot(wi, wh1), Float(m_eta));
            Float T1 = 1.f - R1;

            /* refraction at the first interface */
            vec3 wt = refract(wi, wh1, cos_theta_t1, eta_ti1);
            Float phi_t = dirPhi(wt);
            Float phi_mt = 2.f * phi_t - phi_mi;
            vec3 wmt = sph_dir(-m_tilt, phi_mt);

            /* Simpson's rule weight */
            Float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);

            vec3 wh2;
            Spectrum A_t = exp(sigmaA * 2.f * cos(phi_t - phi_mi) / costheta(wt));
            Float G1 = G(wi, -wt, wmi, wh1);
            if (G1 == 0 || G_(wi, -wt, vec3(wmi.x, 0.f, wmi.z), wh1) == 0)
                continue;

            if (dot(wo, wt) < m_inv_eta - 1e-5f) /* total internal reflection */
                goto TRT;
            wh2 = -wt + m_inv_eta * wo;

            if (dot(wmt, wh2) < 0) /* microfacet invisible from macronormal */
                goto TRT;

            {
                Float rcp_norm_wh2 = rcp((-wt + m_inv_eta * wo).length());
                wh2 = wh2 * rcp_norm_wh2;

                Float dot_wt_wh2 = dot(-wt, wh2);
                Float T2 = 1.f - std::get<0>(fresnel(dot_wt_wh2, Float(m_inv_eta)));
                Float D2 = D(wh2, wmt) * G(-wt, -wo, wmt, wh2);

                /* integrand_of_S_tt / pdf_of_sampling_wt */
                // Spectrum result;
                // if (distr.sample_visible()) {
                //     result = T1 * T2 * smith_g1(-wt, wmi, wh1) * D2 * A_t
                // 	* dot(wi, wmi)
                // 	* dot_wt_wh2 * dot(wo, wh2) * sqr(rcp_norm_wh2)
                // 	* rcp(dot(wt, wmi)) * weight;
                // } else {
                //     result = T1 * T2 * G1 * D2 * A_t * cos_ih1
                // 	* dot_wt_wh2 * dot(wo, wh2) * sqr(rcp_norm_wh2)
                // 	* rcp(dot(wt, wmi)) * weight / dot(wh1, wmi);
                // }
                /* integrand_of_S_tt / pdf_of_sampling_wt */
                bool sampleVisible = true;
                Spectrum result = T1 * T2 * D2 * A_t * dot_wt_wh2 * dot(wo, wh2)
                                  * sqr(rcp_norm_wh2) * rcp(dot(wt, wmi)) * weight *
                                  select(sampleVisible, smith_g1(-wt, wmi, wh1) * dot(wi, wmi),
                                         G1 * cos_ih1 / dot(wh1, wmi));
                if(isinf(luminace(result))) result = Spectrum (0);
             //   masked(result, !enoki::isfinite(result)) = 0;
                S_tt += result;
            }


            TRT:
            vec2 sample2 = m_sampler->getNext2D();
            wh2 = std::get<0>(sample_wh(-wt, wmt, distr, sample2,alphaXY));

            Float cos_th2 = dot(-wt, wh2);
            if (!(cos_th2 > 1e-5f))
                continue;

            /* fresnel coefficient */
            auto [R2, cos_theta_t2, eta_it2, eta_ti2] = fresnel(cos_th2, Float(m_inv_eta));
            vec3 wtr = reflect(wt, wh2);

            Float G2 = G(-wt, -wtr, wmt, wh2);
            if (G2 == 0 || G_(-wt, -wtr, vec3(wmt.x, 0.f, wmt.z), wh2) == 0)
                continue;

            if (dot(-wtr, wo) < m_inv_eta - 1e-5f) /* total internal reflection */
                continue;

            Float phi_tr = dirPhi(wtr);
            Float phi_mtr = phi_mi - 2.f * (phi_t - phi_tr) + M_PI;
            vec3 wmtr = sph_dir(-m_tilt, phi_mtr);

            vec3 wh3 = wtr + m_inv_eta * wo;
            Float G3 = G(wtr, -wo, wmtr, wh3);
            if (dot(wmtr, wh3) < 0 || G3 == 0 || G_(wtr, -wo, vec3(wmtr.x, 0.f, wmtr.z), wh3) == 0)
                continue;

            Float rcp_norm_wh3 = rcp(wh3.length());
            wh3 *= rcp_norm_wh3;

            Float cos_trh3 = dot(wh3, wtr);
            Float T3 = 1.f - std::get<0>(fresnel(cos_trh3, Float(m_inv_eta)));

            Float D3 = D(wh3, wmtr) * G3;
            Spectrum A_tr = exp(sigmaA * 2.f * cos(phi_tr - phi_mt) / costheta(wtr));

            // Spectrum result;
            // if (distr.sample_visible()) {
            // 	result = T1 * R2 * T3 * smith_g1(-wt, wmi, wh1) * smith_g1(-wtr, wmt, wh2)
            // 	    * D3 * cos_trh3 * dot(wi, wmi) * dot(wt, wmt)
            // 	    * dot(wh3, wo)	* sqr(rcp_norm_wh3) * A_t * A_tr * weight
            // 	    / (dot(wt, wmi) * dot(wtr, wmt));
            // } else {
            // 	result = T1 * R2 * T3 * G1 * G2 * D3 * cos_ih1 * cos_th2 * cos_trh3
            // 	    * dot(wh3, -wo)	* sqr(rcp_norm_wh3) * A_t * A_tr * weight
            // 	    / (dot(wh1, wmi) * dot(wh2, wmt) * dot(wt, wmi) * dot(wtr, wmt));
            // }
            bool sample_visible = true;
            Spectrum result = T1 * R2 * T3 * D3 * cos_trh3 * dot(wh3, wo) * sqr(rcp_norm_wh3) *
                              A_t * A_tr * weight / (dot(wt, wmi) * dot(wtr, wmt)) *
                              select(sample_visible,
                                     smith_g1(-wt, wmi, wh1) * smith_g1(-wtr, wmt, wh2) * dot(wi, wmi) *
                                     dot(wt, wmt),
                                     -G1 * G2 * cos_ih1 * cos_th2 / (dot(wh1, wmi) * dot(wh2, wmt)));
            if(isinf(luminace(result))) result = Spectrum (0);
        //    masked(result, !enoki::isfinite(result)) = 0;
            S_trt += result;
        }
        return (S_tt + S_trt) * 1.f / 3.f * res * sqr(m_inv_eta);
    }

#include <fstream>
    Spectrum RoughHair::sampleF(SurfaceEvent & event, const vec2 &sample) const {
        //auto wi = vec3(0, 0.200661987, 0.979660511);
        auto wi  = event.wo;
      //  MTS_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
   //     BSDFSample3f bs = zero<BSDFSample3f>();
   //     Mask active_r, active_tt, active_trt;

        //if (unlikely(!ctx.is_enabled(BSDFFlags::GlossyReflection) || none_or<false>(active)))
          //  return {bs, 0.f};

        /* Construct a microfacet distribution matching the
           roughness values at the current surface position. */
    //    MicrofacetDistribution distr(m_type, m_roughness, m_sample_visible);

        // generate sample
        Float sample_lobe = m_sampler->getNext1D();
        Float sample_h = m_sampler->getNext1D();
        vec2 sample_h1 = sample;
        vec2 sample_h2 = m_sampler->getNext2D();
        vec2 sample_h3 = m_sampler->getNext2D();


        // Float sin_phi_mi = si.dn_du.x;	  /* Use offset h directly from intersection data */
        Float sin_phi_mi = sample_h * 2.f - 1.f;  /* Sample offset h = -sin(phi_m)*/
        Float cos_phi_mi = sqrt(1.f - sqr(sin_phi_mi));
        auto [st, ct] = sincos(m_tilt);
        vec3 wmi(sin_phi_mi * ct, st, cos_phi_mi * ct); /* mesonormal */
        vec3 wmi_(sin_phi_mi, 0.f, cos_phi_mi); /* macronormal */

        if (dot(wmi, wi) < 0 || dot(wmi_, wi) < 0)
            return Spectrum(0); /* macro/mesonormal invisible */

        // sample R lobe
        auto [wh1, pdfh1] = sample_wh(wi, wmi, distr, sample_h1,alphaXY);
        vec3 wr = reflect(wi, wh1);

        bool active_r,active_tt,active_trt;
        bool active = true;
        /* Ensure that this is a valid sample */
        active &= (dot(wr, wh1) > 0 && dot(wr, wmi) > 0 && G_(wi, wr, wmi_, wh1) > 0 && pdfh1 > 0);
        active_r = active;

        /* fresnel coefficient */
        auto [R1, cos_theta_t1, eta_it1, eta_ti1] = fresnel(dot(wi, wh1), Float(m_eta));
        Spectrum R = Spectrum(select(active_r, R1, 0.f));

        // sample TT lobe
        vec3 wt = refract(wi, wh1, cos_theta_t1, eta_ti1);
        Float phi_t = dirPhi(wt);
        Float phi_mi = atan2f(sin_phi_mi, cos_phi_mi);
        Float phi_mt = 2.f * phi_t - phi_mi;
        vec3 wmt = sph_dir(-m_tilt, phi_mt);
        vec3 wmt_ = sph_dir(0, phi_mt);
        auto [wh2, pdfh2] = sample_wh(-wt, wmt, distr, sample_h2,alphaXY);
        vec3 wtr = reflect(wt, wh2);

        /* fresnel coefficient */
        auto [R2, cos_theta_t2, eta_it2, eta_ti2] = fresnel(dot(-wt, wh2), Float(m_inv_eta));

        vec3 wtt = refract(-wt, wh2, cos_theta_t2, eta_ti2);
        active_tt = (active && dot(wt, wh2) < 0 && dot(wmt, wt) < 0 && pdfh2 > 0
                     && G_(-wt, -wtr, vec3(wmt.x, 0.f, wmt.z), wh2) > 0); // visibility
        active_trt = active_tt;
        active_tt &= (dot(wtt, wmt) < 0);
        active_tt &= (cos_theta_t2 != 0); // total internal reflection
        Spectrum T1 = Spectrum (1.f - R1);
        Spectrum T2 = Spectrum (1.f - R2);

        /* absorption */
     //  Spectrum wavelengths = get_spectrum(si);
      //  Spectrum mu_a = fmadd(m_pheomelanin, pheomelanin(wavelengths),
      //                        m_eumelanin * eumelanin(wavelengths));
        Float cos_gamma_t = -cos(phi_t - phi_mi);
        Float cos_theta_wt = sqrt(1.f - sqr(wt.y));
        Spectrum A_t = exp(sigmaA * (-2.f * cos_gamma_t / cos_theta_wt));

        Spectrum TT = select(active_tt, Spectrum(T1 * A_t * T2), Spectrum(0.f));

        // sample TRT lobe
        Float phi_tr = dirPhi(wtr);
        Float phi_mtr = phi_mi - 2.f * (phi_t - phi_tr) + M_PI;
        vec3 wmtr = sph_dir(-m_tilt, phi_mtr);
        vec3 wmtr_ = sph_dir(0, phi_mtr);
        auto [wh3, pdfh3] = sample_wh(wtr, wmtr, distr, sample_h3,alphaXY);

        /* fresnel coefficient */
        auto [R3, cos_theta_t3, eta_it3, eta_ti3] = fresnel(dot(wtr, wh3), Float(m_inv_eta));
        vec3 wtrt = refract(wtr, wh3, cos_theta_t3, eta_ti3);
        active_trt &= (cos_theta_t3 != 0); // total internal reflection
        active_trt &= (dot(wtr, wh3) > 0 && dot(wmtr, wtr) > 0 && dot(wtrt, wmtr) < 0 && pdfh3 > 0
                       && G_(wtr, -wtrt, vec3(wmtr.x, 0.f, wmtr.z), wh3) > 0);
        Spectrum T3 =  Spectrum(1.f - R3);
        Float cos_gamma_t2 = -cos(phi_tr - phi_mt);
        Float cos_theta_wtr = sqrt(1.f - sqr(wtr.y));
        Spectrum A_tr = exp(sigmaA * (-2.f * cos_gamma_t2 / cos_theta_wtr));

        Spectrum TRT = select(active_trt, Spectrum(T1 * R2 * T3 * A_t * A_tr), Spectrum(0));

        // select lobe based on energy
        Float r = average(R);
        Float tt = average(TT);
        Float trt = average(TRT);
        Float total_energy = r + tt + trt;


        active &= (total_energy > 0 && isfinite(total_energy));

        sample_lobe *= total_energy;
        bool selected_r = sample_lobe < r && active_r;
        bool selected_tt = sample_lobe >= r && sample_lobe < (r + tt) && active_tt;
        bool selected_trt = sample_lobe >= (r + tt) && active_trt;

        event.wi = select(selected_r, wr, select(selected_tt, wtt, wtrt));
        if(length(event.wi)>1.1){
            int k = 1;
        }
        //bs = 1.f;
      //  bs.sampled_component = 0;
        event.sampleType = BXDFType(BSDF_GLOSSY | BSDF_REFLECTION);

        Spectrum weight =
                select(selected_r, Spectrum(R / r * total_energy),
                       select(selected_tt, Spectrum(TT / tt * total_energy),
                              select(selected_trt, Spectrum(TRT / trt * total_energy), Spectrum(0.f))));

        Float visibility = select(selected_r, smith_g1(wr, wmi, wh1) * G_(wi, wr, wmi_, wh1),
                                   select(selected_tt, smith_g1(-wt, wmi, wh1) * smith_g1(-wtt, wmt, wh2)
                                                       * G_(wi, -wt, wmi_, wh1) * G_(-wt, -wtt, wmt_, wh2),
                                          select(selected_trt,
                                                 smith_g1(-wt, wmi, wh1) * smith_g1(-wtr, wmt, wh2) *
                                                 smith_g1(-wtrt, wmtr, wh3)
                                                 * G_(wi, -wt, wmi_, wh1) * G_(-wt, -wtr, wmt_, wh2) *
                                                 G_(wtr, -wtrt, wmtr_, wh3),
                                                 0.f)));
        visibility = 1.f;
        // {
        // Float dwh_dwo = select(selected_r, rcp(4.f * dot(wr, wh1)),
        // 		       select(selected_tt,
        // 			      sqr(m_inv_eta) * rcp(squared_norm(-wt + m_inv_eta * wtt)) * dot(-wtt, wh2),
        // 			      select(selected_trt,
        // 				     sqr(m_inv_eta) * rcp(squared_norm(wtr + m_inv_eta * wtrt)) * dot(-wtrt, wh3),
        // 				     0.f)));

        // bs.pdf = abs(dwh_dwo) *
        //     select(selected_r, r / total_energy * pdfh1,
        // 	   select(selected_tt, tt / total_energy * pdfh2,
        // 		  select(selected_trt, trt / total_energy * pdfh3, 0.f)));
        // }

        weight *= visibility;

        /* correction of the cosine foreshortening term */
        // weight *= dot(wi, wmi) / dot(wi, wmi_);

        /* ensure the same pdf is returned for BSDF and emitter sampling */
       return f(event);
    }


    Spectrum RoughHair::f(const SurfaceEvent & event) const {
        auto out = event.wo;
        auto in = event.wi;
   return evalTTAndTRT(out,in) + evalR(out,in);
      //  return rRes + ttAndTrTRes;
    }

RoughHair::RoughHair()  : BSDF(BXDFType(BSDF_REFLECTION | BSDF_GLOSSY)){
    m_sampler = new UniformSampler();
    distr = new GGX();
    m_roughness = 0.135000005;
    m_roughness_squared = sqr(m_roughness);
    m_tilt = 0;
    m_tan_tilt = tan(m_tilt);
    //   m_sampler = new IndependentSampler();
    alphaXY = vec2(m_roughness,m_roughness);
    m_eta =1.54857099;
    m_inv_eta = rcp(m_eta);
    sigmaA = Spectrum(0.308684766, 0.467364371, 0.89055419);
    sigmaA = Spectrum(0.54, 0.52, 1.36);
    m_analytical = true;
}