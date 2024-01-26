#include "Warp.hpp"

/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

vec2 Warp::squareToUniformSquare(const vec2& sample) {
    return sample;
}

Float Warp::squareToUniformSquarePdf(const vec2& sample) {
    return AllInRange<2>(sample, 0.0, 1.0) ? 1.0 : 0.0;
}

static Float TentInverse(Float x) {
    if (x <= .5f)
        return std::sqrt(2 * x) - 1;
    return 1 - std::sqrt(2 - 2 * x);
}

vec2 Warp::squareToTent(const vec2& sample) {
    vec2 res(TentInverse(sample[0]), TentInverse(sample[1]));
    return res;
    //    throw NoriException("Warp::squareToTent() is not yet implemented!");
}

Float Warp::squareToTentPdf(const vec2& p) {
    return (1.f - abs(p[0])) * (1.f - abs(p[1]));
}

vec2 Warp::ConcentricSampleDisk(const vec2& u) {
    vec2 uOffset = 2.f * u - vec2(1, 1);

    // Handle degeneracy at the origin
    if (uOffset.x == 0 && uOffset.y == 0) return vec2(0, 0);

    // Apply concentric mapping to point
    Float theta, r;
    if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
        r     = uOffset.x;
        theta = Constant::PiOver4 * (uOffset.y / uOffset.x);
    } else {
        r     = uOffset.y;
        theta = Constant::PiOver2 - Constant::PiOver4 * (uOffset.x / uOffset.y);
    }
    return r * vec2(std::cos(theta), std::sin(theta));
}

vec2 Warp::squareToUniformDisk(const vec2& sample) {
    auto phi = 2 * sample.x * Constant::PI;
    auto r   = sqrt(sample.y);
    return {r * cos(phi), r * sin(phi)};
}

Float Warp::squareToUniformDiskPdf(const vec2& p) {
    return length(p) < 1.f ? Constant::INV_PI : .0f;
}

vec3 Warp::squareToUniformSphere(const vec2& sample) {
    Float z   = 1 - 2 * sample[0];
    Float r   = std::sqrt(std::fmax((Float)0, (Float)1 - z * z));
    Float phi = 2 * Constant::PI * sample[1];
    return {r * std::cos(phi), r * std::sin(phi), z};
}

Float Warp::squareToUniformSpherePdf(const vec3& v) {
    //    throw NoriException("Warp::squareToUniformSpherePdf() is not yet implemented!");
    return 0.25f * Constant::INV_PI;
}

vec3 Warp::squareToUniformHemisphere(const vec2& sample) {
    Float z   = 1 - 2 * sample[0];
    Float r   = std::sqrt(std::fmax((Float)0, (Float)1 - z * z));
    Float phi = 2 * Constant::PI * sample[1];
    return {r * std::cos(phi), r * std::sin(phi), abs(z)};
    //throw NoriException("Warp::squareToUniformHemisphere() is not yet implemented!");
}

Float Warp::squareToUniformHemispherePdf(const vec3& v) {
    return v[2] >= 0 ? 0.5f * Constant::INV_PI : .0f;
    //    throw NoriException("Warp::squareToUniformHemispherePdf() is not yet implemented!");
}

vec3 Warp::squareToCosineHemisphere(const vec2& sample) {
    Float z   = sqrt(1 - sample.x);
    Float phi = sample.y * 2 * Constant::PI;

    return {sqrt(sample.x) * cos(phi), sqrt(sample.x) * sin(phi), z};
}

Float Warp::squareToCosineHemispherePdf(const vec3& v) {
    return v[2] >= 0 ? v.z * Constant::INV_PI : .0f;
}

//
vec3 Warp::squareToBeckmann(const vec2& sample, Float alpha) {
    auto tan2theta = -alpha * alpha * log(sample.x);
    auto cosTheta  = sqrt(1 / (1 + tan2theta));
    auto sinTheta  = sqrt(1 - cosTheta * cosTheta);
    auto phi       = sample.y * 2 * Constant::PI;
    vec3 t1        = vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
    return t1;
}

//
Float Warp::squareToBeckmannPdf(const vec3& m, Float alpha) {
    if (m.z <= 0)
        return 0.0f;
    auto  cosTheta     = m.z;
    auto  sinTheta     = sqrt(1 - cosTheta * cosTheta);
    auto  tan2Theta    = (sinTheta * sinTheta) / (cosTheta * cosTheta);
    Float azimuthal    = Constant::INV_PI;
    Float longitudinal = exp(-tan2Theta / (alpha * alpha)) / (alpha * alpha * pow(cosTheta, 3));
    return azimuthal * longitudinal;
}
////    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
//}
//    vec3 Warp::squareToBeckmann(const vec2 &sample, Float alpha) {
//        Float phi = 2*Constant::PI*sample[0];
//        Float invAlpha2 = 1.f/alpha; invAlpha2 *= invAlpha2;
//        Float theta = acos(sqrt(1/(1-alpha*alpha*log(sample[1]))));
//
//        return vec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
//    }

//    Float Warp::squareToBeckmannPdf(const vec3 &m, Float alpha) {
//
//
//
//           if(m.z<=0)
//        return 0.0f;
//        if(Frame::cosTheta(m) <= 0)
//            return 0.f;
//            auto cosTheta=m.z;
//            auto sinTheta=sqrt(1-cosTheta*cosTheta);
//            auto tan2Theta=(sinTheta* sinTheta)/(cosTheta*cosTheta);
//          auto t1  =INV_PI ;
//          Float t3=exp(-tan2Theta/(alpha*alpha))  / (alpha*alpha*pow(cosTheta,3));
//            t1*=t3;
//            return t1;
//         cosTheta = Frame::cosTheta(m);
//        Float tanTheta = Frame::tanTheta(m);
//
//        if(cosTheta <= 0)
//            return 0.f;
//
//        Float azimuthal =   INV_PI;
//        auto tan2Theta2=-tanTheta*tanTheta;
//        Float longitudinal = exp(tan2Theta2/(alpha*alpha)) /(alpha*alpha*pow(cosTheta,3));
//
//        auto t2= azimuthal * longitudinal;
//        if(t1!=t2){
//            t2=t1;
//        }
//        return t1;
//    }