#pragma  once

#include "Ray/Intersection.hpp"
#include "Ray/Ray.hpp"
#include "Common/Frame.hpp"
#include "Bsdfs/BsdfTypes.hpp"

struct SurfaceScatterEvent{
    const Intersection *  its;
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

    SurfaceScatterEvent(const SurfaceScatterEvent & event): its(new Intersection(*event.its)),
                                                            frame(event.frame),wo(event.wo),wi(event.wi),sampleType(event.sampleType),
                                                            value(event.value),pdf(event.pdf),flippedFrame(event.flippedFrame)
    {

    }

    SurfaceScatterEvent(){}
    ~SurfaceScatterEvent(){
        delete its;
    }

    Ray sctterRay(const Ray & ray){
        vec3 dir = toWorld(wi);
        vec3 offsetPos = its->p + dir * Constant::EPSILON;
        return Ray(offsetPos,dir,0);
    }
};

