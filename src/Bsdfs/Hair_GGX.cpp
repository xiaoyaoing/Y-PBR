//#include "Hair.h"
//#include "Fresnel.hpp"
//#include <iostream>
//
//std::pair<vec3, Float> sample_wh(const vec3 &wi, const vec3 &wm,
//                                      const MicrofacetDistribution * distr,
//                                      const vec2 &sample1,const vec2 & alphaXY) {
//    /* Coordinate transformation for microfacet sampling */
//    Frame wmFrame;
//    wmFrame.n = wm;
//    wmFrame.tangent = cross(vec3(0.f, 1.f, 0.f), wm);
//    wmFrame.bitTangent = cross(wmFrame.n, wmFrame.tangent);
//    vec3 wh = distr->Sample_wh(wmFrame.toLocal(wi), sample1, alphaXY);
////    if(dot(wh,wmFrame.toLocal(wi))<0)
////        wh = -wh;
//    // wh  = vec3(wh[2],wh[0],wh[1]);
//    auto pdf = distr->Pdf(wmFrame.toLocal(wi), wh, alphaXY);
//    return {wmFrame.toWorld(wh), pdf};
//
//}
//
//vec3 wm_local(const vec3 & wh,const vec3 & wm){
//    Frame wmFrame;
//    wmFrame.n = wm;
//    wmFrame.tangent = cross(vec3(0.f, 1.f, 0.f), wm);
//    wmFrame.bitTangent = cross(wmFrame.n, wmFrame.tangent);
//    return wmFrame.toLocal(wh);
//}
//
//
//
//vec3 wm_world(const vec3 & wh,const vec3 & wm){
//    Frame wmFrame;
//    wmFrame.n = wm;
//    wmFrame.tangent = cross(vec3(0.f, 1.f, 0.f), wm);
//    wmFrame.bitTangent = cross(wmFrame.n, wmFrame.tangent);
//    return wmFrame.toWorld(wh);
//}
//
//std::pair<Float,Float> sincos(Float  angle)  {
//    return {sin(angle), cos(angle)};
//}
//
//inline vec3 sphDir(Float phi,Float theta){
//    auto [sin_theta, cos_theta] = sincos(theta);
//    auto [sin_gamma,   cos_gamma]   = sincos(phi);
//    return vec3 (sin_gamma * cos_theta, sin_theta, cos_gamma * cos_theta);
//}
//
//
//static Float trigInverse(Float x) {
//    return std::min(std::sqrt(std::max(1.0f - x*x, 0.0f)), 1.0f);
//
//}
//
//static Float I0(Float x) {
//    Float result = 1.0f;
//    Float xSq = x * x;
//    Float xi = xSq;
//    Float denom = 4.0f;
//    for ( int i = 1 ; i <= 10 ; ++ i ) {
//        result += xi / denom;
//        xi *= xSq;
//        denom *= 4.0f * Float(( i + 1 ) * ( i + 1 ));
//    }
//    return result;
//}
//
//static Float getH(const SurfaceEvent & event){
//    return std::fabs(event.its->uv.y);
//}
//
//
//static Float Phi(Float gammaI, Float gammaT, int p) {
//    return 2.0f * p * gammaT - 2.0f * gammaI + p * Constant::PI;
//}
//
//
//static Float Phi_1(Float gammaI, Float gammaT, int p) {
//    return p * gammaT -   gammaI + p * Constant::PI;
//}
//
//
//
//static Float Gussian(Float beta, Float theta) {
//    return std::exp(- theta * theta / ( 2.0f * beta * beta )) / ( std::sqrt(2.0f * Constant::PI) * beta );
//}
//
//
//inline Float Logistic(Float x, Float s) {
//    x = std::abs(x);
//    return std::exp(- x / s) / ( s * sqr(1 + std::exp(- x / s)) );
//}
//
//inline Float LogisticCDF(Float x, Float s) {
//    return 1 / ( 1 + std::exp(- x / s) );
//}
//
//
//inline Float TrimmedLogistic(Float x, Float s, Float a, Float b) {
//    while ( x > Constant::PI ) x -= 2 * Constant::PI;
//    while ( x < - Constant::PI ) x += 2 * Constant::PI;
//    return Logistic(x, s) / ( LogisticCDF(b, s) - LogisticCDF(a, s) );
//}
//
//static Float SampleTrimmedLogistic(Float u, Float s, Float a, Float b) {
//    Float k = LogisticCDF(b, s) - LogisticCDF(a, s);
//    Float x = - s * std::log(1 / ( u * k + LogisticCDF(a, s) ) - 1);
//    return clamp(x, a, b);
//}
//
//static std::array < Spectrum, pMax + 1 > Ap(Float cosThetaO, Float eta, Float h,
//                                            const Spectrum & T) {
//    std::array < Spectrum, pMax + 1 > ap;
//    // Compute $p=0$ attenuation at initial cylinder intersection
//    Float cosGammaO = trigInverse(1 - h * h);
//    Float cosTheta = cosThetaO * cosGammaO;
//    Float f = Fresnel::dielectricReflectance(1 / eta, cosTheta);
//    ap[0] = Spectrum(f);
//    // Compute $p=1$ attenuation term
//    ap[1] = sqr(1 - f) * T;
//    // Compute attenuation terms up to $p=_pMax_$
//    for ( int p = 2 ; p < pMax ; ++ p ) ap[p] = ap[p - 1] * T * f;
//    // Compute attenuation term accounting for remaining orders of scattering
//    // ap[pMax] = ap[pMax - 1] * f * T / ( Spectrum(1.f) - T * f );
//    return ap;
//}
//
//static float logI0(float x)
//{
//    if (x > 12.0f)
//        // More stable evaluation of log(I0(x))
//        // See also https://publons.com/discussion/12/
//        return x + 0.5f*(std::log(1.0f/(Constant::TWO_PI*x)) + 1.0f/(8.0f*x));
//    else
//        return std::log(I0(x));
//}
//
//Float Hair::M(Float v, Float sinThetaO, Float sinThetaI, Float cosThetaO, Float cosThetaI) const {
//    Float b = sinThetaO * sinThetaI / v;
//    Float a = cosThetaI * cosThetaO / v;
//    Float csch = 2 / ( exp(1 / v) - exp(- 1 / v) );
//    if(v<0.1)
//        return std::exp(-b + logI0(a) - 1.0f/v + 0.6931f + std::log(1.0f/(2.0f*v)));
//
//    return csch * exp(b) * I0(- a);
//}
//
//Float Hair::NR(Float beta, Float cosO, Float phi, Float h) const {
//    Float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
//    Float deltaPhi = phi + 2.0f * gammaI;
//    deltaPhi = std::fmod(deltaPhi, Constant::TWO_PI);
//    if ( deltaPhi < 0.0f )
//        deltaPhi += Constant::TWO_PI;
//    // return D(beta, deltaPhi) * Fresnel::dielectricReflectance(1.0f / _eta, cosO * cos(gammaI));
////    return TrimmedLogistic(deltaPhi, beta, - Constant::PI, Constant::PI) *
//    return        Fresnel::dielectricReflectance(1.0f / _eta, cosO * cos(gammaI));
//}
//
//vec3 Hair::NP(Float beta, Float cosThetaD, Float phi, int p, Float h) const {
//    Float iorPrime = std::sqrt(_eta * _eta - ( 1.0f - cosThetaD * cosThetaD )) / cosThetaD;
//    Float cosThetaT = std::sqrt(1.0f - ( 1.0f - cosThetaD * cosThetaD ) * sqr(1.0f / _eta));
//    vec3 sigmaAPrime = _sigmaA / cosThetaT;
//    Float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
//    Float gammaT = std::asin(clamp(h / iorPrime, - 1.0f, 1.0f));
//    Float l = 2.0f * std::cos(gammaT);
//
//    Float f = Fresnel::dielectricReflectance(1.0f / _eta, cosThetaD * trigInverse(h));
//    vec3 T = exp(- sigmaAPrime * l);
//    vec3 Aph = ( 1.0f - f ) * ( 1.0f - f ) * T;
//    for ( int i = 1 ; i < p ; ++ i )
//        Aph *= f * T;
//    return Aph;
//    Float deltaPhi = phi - Phi(gammaI, gammaT, p);
//    deltaPhi = std::fmod(deltaPhi, Constant::TWO_PI);
//    if ( deltaPhi < 0.0f )
//        deltaPhi += Constant::TWO_PI;
//    // return Aph * D(beta, deltaPhi);
//    return Aph * TrimmedLogistic(deltaPhi, beta, - Constant::PI, Constant::PI);
//}
//static Float getPhi(vec3 d){
//    return std::atan2(d.x,d.z);
//}
//
//
//
//static Float eval_dielectric(vec3 wi,vec3 wo,vec3 wh,vec3 wm,const GGX & ggx,vec2 alphaXY,Float eta = 1.55){
//    if( dot(wo,wm) < 0)
//    {
//        auto d1 = dot(wi,wm);
//        auto d2= dot(wo,wm);
//        auto phi1 = getPhi(wi);
//        auto phi2 = getPhi(wm);
//        auto phi3= getPhi(wo);
//        return 0;
//    }
//    auto t1 = getPhi(wi);
//    auto t2 = getPhi(wh);
//    auto t3 =getPhi(wo);
//    if(  dot(wh,wi)<0){
//        return 0;
//    }
//    if(dot(wm,wh)<0)
//    {   auto t1= getPhi(wm);
//        auto t2 = getPhi(wh);
//        return 0;}
//    Float whDotIn =  dot(wh,wi);
//    Float whDotOut = dot(wh,wo);
//    Float sqrtDeom = eta * whDotOut  +  whDotIn;
//    return  ggx.D(wm_local(wh,wm),alphaXY) * ggx.G(wm_local(wo,wm), wm_local(wi,wm),alphaXY)  * std::abs(
//            whDotIn * whDotOut  /
//            (dot(wo,wm) * sqrtDeom * sqrtDeom));
//}
//
//static Float pdf_dielectric(vec3 wi,vec3 wo,vec3 wh,vec3 wm,GGX & ggx,vec2 alphaXY,Float eta = 1.55){
//    Float whPdf = ggx.Pdf(wm_local(wo,wm), wm_local(wo,wh),alphaXY);
//    Float sqrtDenom = dot(wo, wh) * eta +  dot(wi, wh);
//    Float dWhDWi =
//            std::abs( dot(wi, wh)) / (sqrtDenom * sqrtDenom);
//    return whPdf * dWhDWi;
//}
//
//template<class T>
//inline  T select(bool mask,T a,T b){
//    return mask?a:b;
//}
//
//std::tuple<Float, Float, Float, Float> fresnel(Float cos_theta_i, Float eta) {
//    bool outside_mask = cos_theta_i > 0;
//    Float  rcp_eta = 1.f/eta,
//            eta_it = select(outside_mask, eta, rcp_eta),
//            eta_ti = select(outside_mask, rcp_eta, eta);
//    Float  cosThetaT;
//    auto r = Fresnel::dielectricReflectance(1/eta,cos_theta_i,cosThetaT);
//    return {r,cosThetaT,eta_it,eta_ti};
//}
//
//static inline vec3 refract(const vec3 &out, vec3 &wh, Float cosThetaT,
//                            Float eta) {
//    auto whDotOut = dot(out,wh);
//    return  (eta * whDotOut - (whDotOut>0?1:-1)*cosThetaT)* wh - eta* out ;
//}
//
//
//
//Spectrum Hair::f(const SurfaceEvent & event) const {
//    Float sinThetaO = event.wo.y;
//    Float thetaO = std::asin(clamp(sinThetaO, - 1.0, 1.0));
//
//    Float sinThetaI = event.wi.y;
//    Float cosThetaI = sqrt(1 - sinThetaI * sinThetaI);
//    Float thetaI = std::asin(clamp(sinThetaI, - 1.0, 1.0));
//
//    Float thetaD = ( thetaI - thetaO ) * 0.5;
//    Float cosThetaD = cos(thetaD);
//
//    Float thetaOR = thetaO - 2.0 * _scaleAngle;
//    Float thetaOTT = thetaO + _scaleAngle;
//    Float thetaOTRT = thetaO + 4 * _scaleAngle;
//
//    Float MR = M(_vR, sin(thetaOR), sinThetaI, cos(thetaOR), cosThetaI);
//    Float MTT = M(_vTT, sin(thetaOTT), sinThetaI, cos(thetaOTT), cosThetaI);
//    Float MTRT = M(_vTRT, sin(thetaOTRT), sinThetaI, cos(thetaOTRT), cosThetaI);
//
//    Float phi = std::atan2(event.wi.x, event.wi.z);
//    if ( phi < 0.0f )
//        phi += Constant::TWO_PI;
//
//    Float h = getH(event);
//    //   h = event.its->uv.y;
//
//    vec3 Nr = vec3(NR(_betaR, trigInverse(event.wo.y), phi, h));
//    vec3 Ntt = NP(_betaR, cosThetaD, phi, 1, h);
//    vec3 Ntrt = NP(_betaR, cosThetaD, phi, 2, h);
//
//    Float iorPrime = std::sqrt(_eta * _eta - ( 1.0f - cosThetaD * cosThetaD )) / cosThetaD;
//    Float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
//    Float gammaT = std::asin(clamp(h / iorPrime, - 1.0f, 1.0f));
//
//    auto wi = event.wi,wo = event.wo;
//
//    auto phiO = std::atan2(wo.x,wo.z);
//    auto wm_r_phi = phiO - gammaI;
//  //  wm_r_phi = phiO  + gammaI - Constant::PI;
//    auto wm_tt_phi = gammaI - (Constant::PI - 2 * gammaT);
//    auto wm_trt_phi = gammaI - 2 *(Constant::PI - 2 * gammaT);
//    wm_tt_phi = Constant::PI + (wm_r_phi - (Constant::PI - 2 * gammaT));
//    wm_trt_phi = Constant::PI+(-wm_tt_phi - (Constant::PI - 2 * gammaT));
//
//
//    auto wm_r_theta =   +   2 * _scaleAngle;
//    auto wm_tt_theta =(  - _scaleAngle);
//    auto wm_trt_theta =(+ 4 * _scaleAngle);
//
//    if(wo.y<0)
//    {
//        wm_r_theta = - wm_r_theta;
//        wm_tt_theta = -wm_tt_theta;
//        wm_trt_theta = - wm_trt_theta;
//    }
//
//    auto wm_r = sphDir(wm_r_phi,wm_r_theta);
//    auto wm_tt = sphDir(wm_tt_phi,wm_tt_theta);
//    auto wm_trt = sphDir(wm_trt_phi,wm_trt_theta);
//
//
//
//
//
//    //auto wo_tt = getSmoothDir(wo,gammaI,gammaT,1);//
// //   auto wo_trt = getSmoothDir(wo,gammaI,gammaT,2);//
//
//    auto [R1, cos_theta_t1, eta_it1, eta_ti1] = fresnel(dot(wo,wm_r), Float(_eta));
//    auto t = getPhi(wm_r);
//    Float whDotOut = dot(wo, wm_r);
//    Float  cosThetaT;
//    auto eta_i = whDotOut>0?1/_eta:_eta;
//    Float  F = Fresnel::dielectricReflectance(1/_eta,whDotOut,cosThetaT);
//    auto wh_r = normalize(wi + wo);
//
//    auto r_res = Spectrum(Fresnel::dielectricReflectance(1/_eta,dot(wh_r,wi)) * ggx.D(wm_local(wh_r,wm_r),alpha_r) * ggx.G(wm_local(wo,wm_r), wm_local(wi,wm_r),alpha_r))/( 4 * abs(dot(wm_r,wo)));
//
//    if(F==1)
//        return r_res;
//    auto wo_tt =  -refract(wo, wm_r, cos_theta_t1, eta_ti1);
//    wo_tt = -((eta_i * whDotOut - (whDotOut>0?1:-1)*cosThetaT)* wm_r - eta_i* wo);
//    auto wo_trt = -Reflect(-wo_tt,wm_tt);
//
//
//
//    auto wh_tt = normalize(wi + wo_tt * _eta);
//    auto wh_trt = normalize(wi + wo_trt * _eta);
//
//    if(dot(wm_tt,wo_tt)>0){
//        int k  =1;
//    }
//
//   // return wi.y * wo.y>0?Spectrum(1,0,0):Spectrum (0,1,0)
//    auto tt_res = Ntt;
//   // return Spectrum(wo.y);
//    return tt_res * Spectrum (eval_dielectric(wi,wo_tt,wh_tt,wm_tt,ggx,alpha_tt));
//        tt_res   *= eval_dielectric(wi,wo_tt,wh_tt,wm_tt,ggx,alpha_tt);
//            auto d  =dot(wm_trt,wo_trt);
//            auto d1 =dot(wm_tt,wo_tt);
//            auto d2 = dot(wm_tt,wi);
//            auto d3= dot(wi,wm_trt);
//   auto phi1 = getPhi(wm_tt);
//   auto phi2 = getPhi(wo_tt);
//   auto phi_r = getPhi(wm_r);
//   auto phi_trt_m = getPhi(wm_trt);
//   auto phi_trt_o = getPhi(wo_trt);
////    if(d1<0)
////        throw("error");
//    //return tt_res;
//
//    auto trt_res = Ntrt * eval_dielectric(wi,wo_trt,wh_trt,wm_trt,ggx,alpha_trt);
//
//
//
////    if(dot(wm_r,wi)<0 || dot(wm_r,wo)<0 || dot(wm_r,wh_r)<0){
////        r_res = Spectrum(0);
////    }
////    auto tt_res = Ntt * eval_dielectric(wi,wo_tt,wh_tt,)
////
////    auto tt_res = Ntt * ggx.D(wm_local(wh_real,wh_tt),alpha_tt) * ggx.G(wm_local(getSmoothDir(wo,gammaI,gammaT,1),wh_r), wm_local(wi,wh_r),alpha_tt);
////    auto trt_res = Ntrt * ggx.D(wm_local(wh_real,wh_trt),alpha_trt) * ggx.G(wm_local(getSmoothDir(wo,gammaI,gammaT,2),wh_r), wm_local(wi,wh_r),alpha_trt);
//    //return trt_res;
//    auto res =  trt_res;
//   return tt_res;
//  //  return Spectrum(0);
//    res+=tt_res + r_res;
//    return res;
////    return   ;
//}
//
//
//Float Hair::Pdf(const SurfaceEvent & event) const {
//
//
//    Float sinThetaO = event.wo.y;
//    Float costhetaO = trigInverse(sinThetaO);
//    Float thetaO = std::asin(clamp(sinThetaO, - 1, 1));
//
//    Float sinThetaI = event.wi.y;
//    Float cosThetaI = sqrt(1 - sinThetaI * sinThetaI);
//    Float thetaI = std::asin(clamp(sinThetaI, - 1.0, 1.0));
//    //First choose a lobe to sample
//    Float h = getH(event);
//    std::array < Float, pMax + 1 > apPdf = ComputeApPdf(costhetaO, h);
//
//    Float thetaOR = thetaO - 2.0 * _scaleAngle;
//    Float thetaOTT = thetaO + _scaleAngle;
//    Float thetaOTRT = thetaO + 4 * _scaleAngle;
//
//    Float phiI = std::atan2(event.wi.x, event.wi.z);
//    Float phiO = std::atan2(event.wo.x, event.wo.z);
//    Float phi = phiI - phiO;
//
//    Float cosThetaD = cos(0.5 * ( std::asin(sinThetaI) - thetaO ));
//    Float iorPrime = std::sqrt(_eta * _eta - ( 1.0f - cosThetaD * cosThetaD )) / cosThetaD;
//    Float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
//    Float gammaT = std::asin(clamp(h / iorPrime, - 1.0f, 1.0f));
//
//    Float pdf = 0;
//
//
//    auto wm_r_phi =  gammaI;
//    auto wm_tt_phi = gammaI - (Constant::PI - 2 * gammaT);
//    auto wm_trt_phi = gammaI - 2 *(Constant::PI - 2 * gammaT);
//
//    auto wm_r_theta = 2 * _scaleAngle;
//    auto wm_tt_theta = - _scaleAngle;
//    auto wm_trt_theta = 4 * _scaleAngle;
//
//    auto wh_r = sphDir(wm_r_phi,wm_r_theta);
//    auto wh_tt = sphDir(wm_tt_phi,wm_tt_theta);
//    auto wh_trt = sphDir(wm_trt_phi,wm_trt_theta);
//
//    auto wh_real = normalize((event.wo + event.wi));
//    auto wo = event.wo;
//    auto wi = event.wi;
//    auto pdf_r = apPdf[0] * ggx.Pdf(wm_local(wh_r,wo), wm_local(wh_r,wi),alpha_r)/(4 * dot(wo,wh_r));
//    auto pdf_tt = apPdf[1] * ggx.Pdf(wm_local(wh_tt,getSmoothDir(wo,gammaI,gammaT,1)), wm_local(wh_tt,wi),alpha_tt)/(4 * dot(getSmoothDir(wo,gammaI,gammaT,1),wh_tt));
//    auto pdf_trt= apPdf[2] * ggx.Pdf(wm_local(wh_trt,getSmoothDir(wo,gammaI,gammaT,2)), wm_local(wh_trt,wi),alpha_trt)/(4 * dot(getSmoothDir(wo,gammaI,gammaT,2),wh_trt));
//
//    return pdf_r + pdf_tt + pdf_trt;
//
//
//    pdf += M(_vR, sin(thetaOR), sinThetaI, cos(thetaOR), cosThetaI) * apPdf[0] * TrimmedLogistic(phi- Phi(gammaI,gammaT,0),_betaR,-Constant::PI,Constant::PI);
//    pdf += M(_vTT, sin(thetaOTT), sinThetaI, cos(thetaOTT), cosThetaI) * apPdf[1] * TrimmedLogistic(phi- Phi(gammaI,gammaT,1),_betaR,-Constant::PI,Constant::PI);
//    pdf += M(_vTRT, sin(thetaOTRT), sinThetaI, cos(thetaOTRT), cosThetaI) * apPdf[2] *TrimmedLogistic(phi- Phi(gammaI,gammaT,2),_betaR,-Constant::PI,Constant::PI);
//    return pdf;
//}
//
//
//static inline vec3  reflect(const vec3 & out,const vec3 & wh){
//    return -out + 2 * dot(wh,out) * wh;
//}
//
//Spectrum Hair::sampleF(SurfaceEvent & event, const vec2 & u) const {
//    //Hair-samplineg requires 4 randoms.
//    vec2 u0 = DemuxFloat(u[0]), u1 = DemuxFloat(u[1]);
//    Float sinThetaO = event.wo.y;
//    Float costhetaO = trigInverse(sinThetaO);
//    Float thetaO = std::asin(clamp(sinThetaO, - 1, 1));
//    //First choose a lobe to sample
//    Float h = getH(event);
//    std::array < Float, pMax + 1 > apPdf = ComputeApPdf(costhetaO, h);
//    int p;
//    for ( p = 0 ; p < pMax ; ++ p ) {
//        if ( u0[0] <= apPdf[p] )
//            break;
//        u0[0] -= apPdf[p];
//    }
//    Float theta, v;
//
//
//    Float phiO = std::atan2(event.wo.x, event.wo.z);
//    Float deltaphi;
//    Float iorPrime = _eta;
//    Float gammaI = std::asin(clamp(h, - 1.0f, 1.0f));
//    Float gammaT = std::asin(clamp(h / iorPrime, - 1.0f, 1.0f));
//
//
//
//    auto wm_r_phi =  gammaI;
//    auto wm_tt_phi = gammaI - (Constant::PI - 2 * gammaT);
//    auto wm_trt_phi = gammaI - 2 *(Constant::PI - 2 * gammaT);
//
//    auto wm_r_theta = 2 * _scaleAngle;
//    auto wm_tt_theta = - _scaleAngle;
//    auto wm_trt_theta = 4 * _scaleAngle;
//
//    auto wm_r = sphDir(wm_r_phi,wm_r_theta);
//    auto wm_tt = sphDir(wm_tt_phi,wm_tt_theta);
//    auto wm_trt = sphDir(wm_trt_phi,wm_trt_theta);
//
//    vec3 wi,wo = event.wo;
//    auto wo_0 = getSmoothDir(wo,gammaI,gammaT,0);
//    auto wo_1 = getSmoothDir(wo,gammaI,gammaT,1);
//    auto wo_2 = getSmoothDir(wo,gammaI,gammaT,2);
//
//
//    if ( p == 0 ) {
//        auto wh_r = sample_wh(wo_0,wm_r,&ggx,u1,alpha_r).first;
//        wi = reflect(wo,wh_r);
//    } else if ( p == 1 ) {
//       // wi = reflect(getSmoothDir(wo,gammaI,gammaT,1),wh_tt);
//    } else if ( p == 2 ) {
////        wo = reflect()
//    } else { throw ( "Invalid P" ); }
//
//
//
//
//
//
//    deltaphi = Phi(gammaI, gammaT, p) + SampleTrimmedLogistic(u0[1], _betaR, - Constant::PI, Constant::PI);
//    Float phi = phiO + deltaphi;
//    Float sinPhi = sin(phi), cosPhi = cos(phi);
//
////    event.wi = vec3(sinPhi * cosThetaI, sinThetaI, cosPhi * cosThetaI);
//    event.sampleType = this->m_type;
//    event.pdf = Pdf(event);
//    return f(event);
//}
//
//
//
//
//Float Hair::D(Float beta, Float phi) const {
//    Float result = 0.0f;
//    Float delta;
//    Float shift = 0.0f;
//    do {
//        delta = Gussian(beta, phi + shift) + Gussian(beta, phi - shift - Constant::TWO_PI);
//        result += delta;
//        shift += Constant::TWO_PI;
//    } while ( delta > 1e-4f );
//    return result;
//}
//
//inline Float RoughnessToAlpha(Float roughness) {
//    roughness = std::max(roughness, (Float)1e-3);
//    Float x = std::log(roughness);
//    return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x +
//           0.000640711f * x * x * x * x;
//}
//
//Hair::Hair(const Json & json) : BSDF(BXDFType(BSDF_GLOSSY | BSDF_TRANSMISSION | BSDF_REFLECTION)) {
//    _scaleAngle = getOptional(json, "scale_angle", 2.5);
//    _roughness = getOptional(json, "roughness", 0.3);
//    Float melaninRatio = getOptional(json, "melanin_ratio", 1.f);
//    Float melaninConcentration = getOptional(json, "melanin_concentration", 1.3);
//    bool overrideSigmaA = containsAndGet(json, "sigma_a", _sigmaA);
//    if ( ! overrideSigmaA ) {
//        const vec3 eumelaninSigmaA = vec3(0.419f, 0.697f, 1.37f);
//        const vec3 pheomelaninSigmaA = vec3(0.187f, 0.4f, 1.05f);
//
//        _sigmaA = melaninConcentration * lerp(eumelaninSigmaA, pheomelaninSigmaA, melaninRatio);
//    }
//
//    betaM = getOptional(json,"beta_m",0.3);
//    betaN = getOptional(json,"beta_n",0.3);
//
//    _betaR = std::max(Constant::PI * 0.5f * betaM, 0.04f);
//    _betaTT = _betaR * 0.5f;
//    _betaTRT = _betaR * 2.0f;
//
//    _vR = _betaR * _betaR;
//    _vTT = _betaTT * _betaTT;
//    _vTRT = _betaTRT * _betaTRT;
//
//    _betaR = std::max(Constant::PI * 0.5f * betaN, 0.04f);
//
//
//    _scaleAngle = Angle::degToRad(_scaleAngle);
//    alpha_r = vec2(RoughnessToAlpha(betaM), RoughnessToAlpha(betaN));
//    alpha_tt = alpha_r * 0.5f;
//    alpha_trt = alpha_r * 2.0f;
//
//
//
//}
//
////return sinTheta
//Float Hair::sampleM(Float v, Float sinThetaO, Float cosThetaO, Float xi1, Float xi2) const {
//    Float cosTheta = 1.0f + v * std::log(xi1 + ( 1.0f - xi1 ) * std::exp(- 2.0f / v));
//    Float sinTheta = trigInverse(cosTheta);
//    Float cosPhi = std::cos(Constant::TWO_PI * xi2);
//    return clamp(-cosTheta * sinThetaO + sinTheta * cosPhi * cosThetaO,-1,1);
//}
//
//std::array < Float, pMax + 1 > Hair::ComputeApPdf(Float cosThetaO, Float h) const {
//
//    Float sinThetaO = trigInverse(cosThetaO);
//
//    Float sinThetaT = sinThetaO / _eta;
//    Float cosThetaT = trigInverse(sinThetaT);
//
//    Float etap = std::sqrt(_eta * _eta - sqr(sinThetaO)) / cosThetaO;
//    Float sinGammaT = h / etap;
//    Float cosGammaT = trigInverse(sinGammaT);
//
//    Spectrum T = exp(- _sigmaA * 2.f * cosGammaT / cosThetaT);
//    std::array < Spectrum, pMax + 1 > ap = Ap(cosThetaO, _eta, h, T);
//
//    // Compute $A_p$ PDF from individual $A_p$ terms
//    std::array < Float, pMax + 1 > apPdf;
//    Float sumY =
//            std::accumulate(ap.begin(), ap.end()-1, Float(0),
//                            [](Float s, const Spectrum & ap) { return s + luminace(ap); });
//    for ( int i = 0 ; i < pMax ; ++ i ) apPdf[i] = luminace(ap[i]) / sumY;
//    if( isnan(apPdf[0])){
//
//    }
//    return apPdf;
//}
//
//
////入射方向
//vec3 Hair::getSmoothDir(vec3 wo, Float gammaI, Float gammaT, int p) const{
//    auto theta = asin(wo.y) ;
//    auto phi = std::atan2(wo.x,wo.z);
//    if(p==0)
//        return sphDir(theta-2*_scaleAngle,phi + Phi_1(gammaI,gammaT,p));
//    if(p==1)
//        return sphDir(theta+_scaleAngle,phi + Phi_1(gammaI,gammaT,p));
//    if(p==2)
//        return sphDir(theta+4*_scaleAngle,phi+ Phi_1(gammaI,gammaT,p));
//
//}
//