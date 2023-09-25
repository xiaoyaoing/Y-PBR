#pragma  once

#include <iostream>
#include "Ray/Intersection.hpp"
#include "Ray/Ray.hpp"
#include "Common/Frame.hpp"
#include "Bsdfs/BsdfTypes.hpp"
#include "Mediums/PhaseFunction.hpp"

struct SurfaceEvent {
    const Intersection * its = nullptr;
    Frame frame;
    vec3 wo, wi;
    BXDFType sampleType;
    BXDFType requestType;
    Float value;
    Float pdf;
    bool flippedFrame;
    bool itsCopyied = false;
public:
    vec3 toLocal(const vec3 & w) const {
        return frame.toLocal(w);
    }
    vec3 toWorld(const vec3 & w) const {
        return frame.toWorld(w);
    }
//    SurfaceEvent(const SurfaceEvent & event) : its(event.its),
//                                               frame(event.frame), wo(event.wo), wi(event.wi),
//                                               sampleType(event.sampleType),
//                                               value(event.value), pdf(event.pdf),
//                                               requestType(event.requestType),
//                                               flippedFrame(event.flippedFrame),
//                                               itsCopyied(true){
//    }
    SurfaceEvent() {}

    Ray sctterRay(const vec3 & w) ;
    Ray sctterRay() {
        return sctterRay(toWorld(wi));
    }

    SurfaceEvent makeFlipQuery() const;
    SurfaceEvent makeWarpQuery(const vec3 &newWi, const vec3 &newWo) const;
};

struct VolumeEvent {
    bool exited;
    vec3 rayDir,p;
    const PhaseFunction * phase;
    Float pdf;
};


