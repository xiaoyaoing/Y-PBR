#include "Integrator.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"
class AbstractPathTracer : public SamplerIntegrator{
protected:
    SurfaceScatterEvent makeLocalScatterEvent(const Intersection * its) const  ;
//    AbstractPathTracer(Json j) :Integrator(j){ }
};