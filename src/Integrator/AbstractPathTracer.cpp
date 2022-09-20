#include "AbstractPathTracer.hpp"
#include "Bsdfs/Reflection.hpp"

SurfaceScatterEvent AbstractPathTracer::makeLocalScatterEvent(const Intersection * its) const  {
    SurfaceScatterEvent event;
    Frame frame = its->primitive->setTangentFrame(its);

    bool enableTwoSideShading = true ;
    if(abs(its->w.x +0.027210)<0.01){
        int k=1;
    }
    bool hitBackSide = dot(its->w,its->Ng)>0;
    bool isTransmissive  = its->bsdf->MatchesFlags(BSDF_TRANSMISSION);
    bool flippedFrame = false;
    //todo add this to config class
    if(enableTwoSideShading && hitBackSide && !isTransmissive ){
        frame.n = -frame.n;
        frame.s = -frame.s;
        flippedFrame = true;
    }

    event.frame = frame;
    event.its= its;
    event.wo =event.toLocal(-its->w);
    event.flippedFrame= flippedFrame;

    return event;


}

