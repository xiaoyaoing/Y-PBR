#include "SurfaceScatterEvent.hpp"

Ray SurfaceEvent::sctterRay(const vec3& w) {
    vec3 offsetPos = its->p + w * its->epsilon;
    return Ray(offsetPos, w, 0);
}

SurfaceEvent SurfaceEvent::makeWarpQuery(const vec3& newWi, const vec3& newWo) const {

    SurfaceEvent event(*this);
    event.wo = newWo;
    event.wi = newWi;
    return event;
}

SurfaceEvent SurfaceEvent::makeFlipQuery() const {
    SurfaceEvent event(*this);
    event.wi = wo;
    event.wo = wi;
    return event;
}