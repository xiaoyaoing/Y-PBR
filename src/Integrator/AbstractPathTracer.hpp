#include "Integrator.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"
class AbstractPathTracer : public  Integrator{
protected:
    SurfaceScatterEvent makeLocalScatterEvent(const Intersection * its) const  ;
    AbstractPathTracer(nlohmann::json j) :Integrator(j){ }
};