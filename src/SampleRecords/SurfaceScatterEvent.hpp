#pragma  once

#include "Ray/Intersection.hpp"
#include "Ray/Ray.hpp"
#include "Common/Frame.hpp"
#include "Bsdfs/BsdfTypes.hpp"

struct SurfaceEvent {
    const Intersection * its;
    Frame frame;
    vec3 wo, wi;
    BXDFType sampleType;
    BXDFType requestType;
    Float value;
    Float pdf;
    bool flippedFrame;
public:
    vec3 toLocal(const vec3 & w) const {
        return frame.toLocal(w);
    }
    vec3 toWorld(const vec3 & w) const {
        return frame.toWorld(w);
    }
    SurfaceEvent(const SurfaceEvent & event) : its(new Intersection(* event.its)),
                                               frame(event.frame), wo(event.wo), wi(event.wi),
                                               sampleType(event.sampleType),
                                               value(event.value), pdf(event.pdf),
                                               flippedFrame(event.flippedFrame) {}
    SurfaceEvent() {}
    ~SurfaceEvent( ) {
        //delete its;
    }
    Ray sctterRay(const vec3 & w) {
        vec3 offsetPos = its->p + w * its->epsilon;
        return Ray(offsetPos, w, 0);
    }
};

struct VolumeEvent {
    bool exited;
};

