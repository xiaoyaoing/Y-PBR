#include "SkyDome.hpp"
#include "Colors/Spectral.hpp"
#include "Common/Transform.hpp"
#include "scene.hpp"
#include "Texture/BitMapTexture.hpp"
#include "SampleRecords/PositionAndDirectionSample.h"
#include <skylight/ArHosekSkyModel.h>

Spectrum
SkyDome::sampleLi(const vec3& ref, const vec2& u, vec3* wi, Float* pdf, Float* distance) const {
    Float mapPdf;
    vec2  uv = _sky->sample(MAP_SPHERICAL, u, &mapPdf);
    if (!mapPdf)
        return Spectrum(0);
    float sinTheta;
    *wi = uvToDirection(uv, sinTheta);
    if (sinTheta == 0)
        *pdf = 0;
    else
        *pdf = mapPdf / (2 * Constant::PI * Constant::PI * sinTheta);
    *distance = 2 * _worldRadius;
    return _sky->eval(uv);
}

PositionAndDirectionSample SkyDome::sampleDirect(const vec2& positionSample, const vec2& dirSample) const {
    _NOT_IMPLEMENT_ERROR;
}

//Spectrum SkyDome::directLighting(const Intersection & intr) const {
//    return _sky->eval(directionToUV(intr.w));
//}

Spectrum SkyDome::Le(const Ray& ray) const {
    return _sky->eval(directionToUV(ray.d));
}

Spectrum SkyDome::Power() {
    return 4 * Constant::PI * _worldRadius * _worldRadius * _sky->average();
}

static const int SizeX      = 512;
static const int SizeY      = 256;
static const int NumSamples = 10;

static void fillImage(ArHosekSkyModelState* state, float* lambdas, vec3* weights, vec3* img, vec3 sun, float gammaScale) {
    for (int y = 0; y < SizeY / 2; ++y) {
        float theta = (y + 0.5f) * Constant::PI / SizeY;
        for (int x = 0; x < SizeX; ++x) {
            float phi   = (x + 0.5f) * Constant::TWO_PI / SizeX;
            vec3  v     = vec3(std::cos(phi) * std::sin(theta), std::cos(theta), std::sin(phi) * std::sin(theta));
            float gamma = clamp(std::acos(clamp(dot(v, sun), -1.0f, 1.0f)) * gammaScale, 0.0f, Constant::PI);

            vec3 xyz(0.0f);
            for (int i = 0; i < NumSamples; ++i)
                xyz += weights[i] * Float(arhosekskymodel_radiance(state, theta, gamma, lambdas[i]));

            img[x + y * SizeX] += Spectral::xyzToRgb(xyz);
        }
    }
}

void SkyDome::Preprocess(const Scene& scene) {
    float lambdas[NumSamples];
    vec3  weights[NumSamples];

    Spectral::spectralXyzWeights(NumSamples, lambdas, weights);

    vec3 sun = transformVector(_transform, vec3(0, 1, 0));

    float sunElevation = std::asin(clamp(sun.y, -1.0f, 1.0f));

    ArHosekSkyModelState* sunState = arhosekskymodelstate_alienworld_alloc_init(sunElevation, _intensity, _temperature, _turbidity, 0.2f);

    std::unique_ptr<vec3[]> img(new vec3[SizeX * SizeY]);
    std::memset(img.get(), 0, SizeX * SizeY * sizeof(img[0]));

    fillImage(sunState, lambdas, weights, img.get(), sun, 1.0f);

    for (int y = SizeY / 2; y < std::min(SizeY / 2 + 2, SizeY); ++y)
        std::memcpy(img.get() + y * SizeX, img.get() + (SizeY / 2 - 1) * SizeX, SizeX * sizeof(img[0]));

    _sky = std::make_shared<BitMapTexture<Spectrum>>(img.release(), SizeX, SizeY, BitMapTexture<Spectrum>::TexelType::RGB_HDR, true, false);
    _sky->makeSamplable(MAP_SPHERICAL);

    scene.getWorldBound().BoundingSphere(&_worldCenter, &_worldRadius);
}

Float SkyDome::PdfLi(const Intersection& pShape, const vec3& ref) const {
    Float sinTheta;
    vec2  uv = directionToUV(pShape.w, sinTheta);
    if (sinTheta == 0) return 0;
    return Constant::INV_PI * Constant::INV_TWO_PI * _sky->pdf(MAP_SPHERICAL, uv) / sinTheta;
}

bool SkyDome::isDeltaLight() const {
    return false;
}

vec2 SkyDome::directionToUV(const vec3& wi) const {

    return vec2(std::atan2(wi.z, wi.x) * Constant::INV_TWO_PI + 0.5f, std::acos(wi.y) * Constant::INV_PI);
}

vec2 SkyDome::directionToUV(const vec3& wi, float& sinTheta) const {
    sinTheta = std::sqrt(std::max(1.0f - wi.y * wi.y, 0.0f));
    return vec2(std::atan2(wi.z, wi.x) * Constant::INV_TWO_PI + 0.5f, std::acos(wi.y) * Constant::INV_PI);
}

vec3 SkyDome::uvToDirection(vec2 uv, float& sinTheta) const {
    float phi   = (uv.x - 0.5f) * Constant::TWO_PI;
    float theta = uv.y * Constant::PI;
    sinTheta    = std::sin(theta);
    return vec3(
        std::cos(phi) * sinTheta,
        -std::cos(theta),
        std::sin(phi) * sinTheta);
}

SkyDome::SkyDome(const Json& json) {
    _transform   = getOptional(json, "transform", getIndentifyTransform());
    _temperature = getOptional(json, "temoerature", 5777);
    _gammaScale  = getOptional(json, "gamma_scale", 1);
    _turbidity   = getOptional(json, "turbidity", 3);
    _intensity   = getOptional(json, "intensity", 2);
}