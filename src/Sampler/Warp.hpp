/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include "../Common/math.hpp"

/// A collection of useful warping functions for importance sampling
/// I learned it from Nori  http://wjakob.github.io/nori/#
class Warp {
public:
    /// Dummy warping function: takes uniformly distributed points in a square and just returns them
    static vec2 squareToUniformSquare(const vec2& sample);

    /// Probability density of \ref squareToUniformSquare()
    static Float squareToUniformSquarePdf(const vec2& p);

    /// sample a 2D tent distribution
    static vec2 squareToTent(const vec2& sample);

    /// Probability density of \ref squareToTent()
    static Float squareToTentPdf(const vec2& p);

    /// Uniformly sample a vector on a 2D disk with radius 1, centered around the origin
    static vec2 squareToUniformDisk(const vec2& sample);

    /// Probability density of \ref squareToUniformDisk()
    static Float squareToUniformDiskPdf(const vec2& p);

    static vec2 ConcentricSampleDisk(const vec2& u);

    /// Uniformly sample a vector on the unit sphere with respect to solid angles
    static vec3 squareToUniformSphere(const vec2& sample);

    /// Probability density of \ref squareToUniformSphere()
    static Float squareToUniformSpherePdf(const vec3& v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to solid angles
    static vec3 squareToUniformHemisphere(const vec2& sample);

    /// Probability density of \ref squareToUniformHemisphere()
    static Float squareToUniformHemispherePdf(const vec3& v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to projected solid angles
    static vec3 squareToCosineHemisphere(const vec2& sample);

    /// Probability density of \ref squareToCosineHemisphere()
    static Float squareToCosineHemispherePdf(const vec3& v);

    /// Warp a uniformly distributed square sample to a Beckmann distribution * cosine for the given 'alpha' parameter
    static vec3 squareToBeckmann(const vec2& sample, Float alpha);

    /// Probability density of \ref squareToBeckmann()
    static Float squareToBeckmannPdf(const vec3& m, Float alpha);
};