#include "Primitive.hpp"
#include "scene.hpp"

Intersection Primitive::sample(const vec3& ref, const vec2& u, Float* pdf, vec3* wi) const {
    auto intr = sample(u, pdf);

    *wi = intr.p - ref;
    if (length2(*wi) == 0)
        *pdf = 0;
    else {
        auto l = length(*wi);
        *wi /= l;
        // Convert from area measure, as returned by the Sample() call
        // above, to solid angle measure.
        *pdf *= l * l / abs(dot(*wi, intr.Ng));
        if (std::isinf(*pdf)) *pdf = 0.f;
    }
    return intr;
}

Float Primitive::directPdf(const Intersection& pShape, vec3 ref) const {
    return distance2(pShape.p, ref) * InvArea() / (absDot(pShape.Ng, -pShape.w));
}

bool Primitive::initRTC() {
    _geom = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_USER);

    rtcSetGeometryUserPrimitiveCount(_geom, 1);
    rtcSetGeometryUserData(_geom, this);
    rtcSetGeometryBoundsFunction(_geom, &EmbreeUtils::instanceBoundsFunc, nullptr);
    rtcSetGeometryIntersectFunction(_geom, &EmbreeUtils::instanceIntersectFunc);
    rtcSetGeometryOccludedFunction(_geom, &EmbreeUtils::instanceOccludedFunc);
    rtcCommitGeometry(_geom);

    return true;
}

void Primitive::load(const Json& json, const Scene& scene) {
    bsdf   = scene.fetchBSDF(getOptional(json, "bsdf", Json()));
    bssrdf = scene.fetchBSSRDF(getOptional(json, "bssrdf", std::string()));
    if (json.contains("int_medium")) inMedium = scene.fetchMedium(json["int_medium"]);
    if (json.contains("ext_medium")) outMedium = scene.fetchMedium(json["ext_medium"]);
}