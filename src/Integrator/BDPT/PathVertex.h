#pragma once

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
        VertexRecord(){}
        VertexRecord(const SurfaceRecord & surfaceRecord){

        }
        ~VertexRecord(){}
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
//    PathVertex(PathVertex & other) = delete;
    PathVertex(const Light * light,Spectrum Le,Float pdf){
        type = VertexType::Light;
        _sampler.light = light;
        beta = Le;
        pdfFwd = pdf;
    }
    void initLight(const Light * light,Float lightPdf){
        type = VertexType::Light;
        _sampler.light = light;
        _record.lightRecord.lightPdf = lightPdf;
    }

    PathVertex(const Light * light,Float lightPdf){
        initLight(light,lightPdf);
    }
    void initCamera(const Camera * camera,vec2 point){
        type = VertexType::Camera;
        _sampler.camera = camera;
        _record.cameraRecord.pixel = point;
    }

    PathVertex(const Camera * camera,vec2 point){
       initCamera(camera,point);
    }

    PathVertex(const Camera * camera,PositionAndDirectionSample sample){
        type = VertexType::Camera;
        _sampler.camera = camera;
        _record.cameraRecord.sample = sample;
        beta = sample.weight / (sample.dirPdf * sample.posPdf);
    }

    PathVertex(const Light * light,PositionAndDirectionSample sample,Float lightPdf){
        type = VertexType::Light;
        _sampler.light = light; ;
        _record.lightRecord.sample = sample;
        _record.lightRecord.lightPdf = lightPdf;
        beta = sample.weight/lightPdf;
        pdfFwd =  0;
    }

//    void initSurface(const SurfaceRecord & record,const Spectrum  &beta){
//        type = VertexType::Surface;
//        _record = record;
//        _sampler.bsdf = record.its.bsdf;
//        this->beta = beta;
//    }

    PathVertex(const SurfaceRecord & record,const Spectrum  &beta):type(VertexType::Surface),
    _record(record)
    {
        _sampler.bsdf = record.its.bsdf;
        //avoid pointer error
        //_record.surfaceRecord.event.its = &_record.surfaceRecord.its;
        this->beta = beta;
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
        return false;
    }
};