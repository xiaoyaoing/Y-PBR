#pragma  once

#include "Ray/Intersection.hpp"
#include "Common/Frame.hpp"
#include "Bsdfs/BsdfTypes.hpp"

struct SurfaceScatterEvent{
    const Intersection * its;
    Frame frame;
    vec3 wo,wi;
    BXDFType sampleType;
    Float value;
    Float pdf;
    bool flippedFrame;

public:
    vec3  toLocal(const vec3 & w) const {
        return frame.toLocal(w);
    }

    vec3  toWorld(const vec3 & w) const {
       return frame.toWorld(w);
    }
};

