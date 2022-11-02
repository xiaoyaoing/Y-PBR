#include "Primitive.hpp"


Intersection Primitive::sample(const Intersection & ref, const vec2 & u, Float * pdf) const {
    auto intr = sample(u, pdf);

    vec3 wi = intr.p - ref.p;
    if (length2(wi) == 0)
        *pdf = 0;
    else {
        wi = normalize(wi);
        // Convert from area measure, as returned by the Sample() call
        // above, to solid angle measure.
        *pdf *= length2(ref.p-intr.p) / abs(dot(intr.Ng, -wi));
        if (std::isinf(*pdf)) *pdf = 0.f;
    }
    return intr;
}

RTCGeometry Primitive::initRTC( ) {
    _geom = rtcNewGeometry(EmbreeUtils::getDevice(), RTC_GEOMETRY_TYPE_USER);

    rtcSetGeometryUserPrimitiveCount(_geom, 1);
    rtcSetGeometryUserData(_geom, this);
    rtcSetGeometryBoundsFunction(_geom, &EmbreeUtils::instanceBoundsFunc, nullptr);
    rtcSetGeometryIntersectFunction(_geom, &EmbreeUtils::instanceIntersectFunc);
    rtcSetGeometryOccludedFunction(_geom, &EmbreeUtils::instanceOccludedFunc);
    rtcCommitGeometry(_geom);

}


