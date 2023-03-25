#include "Ray/Intersection.hpp"
#include "Records.h"
#include "SampleRecords/SurfaceScatterEvent.hpp"
enum class VertexType{
    Camera,Light,Surface,Medium
};
class Camera;
class Light;
class Medium;
struct EndPointInter : Intersection{
    union{
        const Camera * camera;
        const Light * light;
    };
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
    VertexRecord record;
    VertexSampler sampler;
public:
    PathVertex(const Light * light,Spectrum Le,Float pdf){
        type = VertexType::Light;
        sampler.light = light;
        beta = Le;
        pdfFwd = pdf;
    }
    PathVertex(const Camera * camera,vec2 point){
        type = VertexType::Camera;
        sampler.camera = camera;
        record.cameraRecord.pos = point;
    }

    PathVertex(const SurfaceEvent & event,const Spectrum  &beta,Float pdf,PathVertex & prev){

    }
};