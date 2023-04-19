#include "Ray/Intersection.hpp"
#include "Records.h"
#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "scene.hpp"
#include "Integrator/TraceHelper.h"
enum class VertexType{
    Camera,Light,Surface,Medium
};
class Camera;
class Light;
class Medium;


struct PathState{
    PathState(Sampler &sampler):sampler(sampler){
        bounce =  0;
    }

    Sampler & sampler;
    int bounce;
    Medium * medium = nullptr;
};

class PathVertex {
public:

    VertexType type;
    Spectrum  beta;
    Float pdfBack,pdfFwd;
    union VertexRecord{
        LightRecord lightRecord;
        CameraRecord cameraRecord;
        SurfaceRecord surfaceRecord;
        MediumRecord mediumRecord;
    };
    union VertexSampler{
        const Light * light;
        const Camera * camera;
        const BSDF * bsdf;
        const Medium * medium;
    };
    VertexRecord _record;
    VertexSampler _sampler;
public:
    PathVertex(){}
    PathVertex(const Light * light,Spectrum Le,Float pdf){
        type = VertexType::Light;
        _sampler.light = light;
        beta = Le;
        pdfFwd = pdf;
    }
    PathVertex(const Light * light,Float lightPdf){
        type = VertexType::Light;
        _sampler.light = light;
        _record.lightRecord.lightPdf = lightPdf;
    }
    PathVertex(const Camera * camera,vec2 point){
        type = VertexType::Camera;
        _sampler.camera = camera;
        _record.cameraRecord.pixel = point;
    }

    PathVertex(const Camera * camera,PositionAndDirectionSample sample){
        type = VertexType::Camera;
        _sampler.camera = camera;
        _record.cameraRecord.sample = sample;
        beta = sample.weight;
    }

    PathVertex(const Light * light,PositionAndDirectionSample sample){
        type = VertexType::Light;
        _sampler.light = light; ;
        _record.cameraRecord.sample = sample;
        beta = sample.weight;
    }

    PathVertex(const SurfaceEvent & event,const Spectrum  &beta):type(VertexType::Surface){
        _sampler.bsdf = event.its->bsdf;
        _record.surfaceRecord.event = event;
    }

    bool sampleNext(const Scene & scene, bool adjoint, PathState &state, PathVertex *prev, PathVertex & next);
    bool sampleRootVertex(PathState &state);

    inline bool isLight() const {
        return type == VertexType::Light;
    }

    bool canConnect() const;

    vec3 pos() const;

    Spectrum  eval(const PathVertex & vertex,bool adjoint) const;

    bool isInfiniteLight() const {

    }
};