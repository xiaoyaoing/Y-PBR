#pragma once

#include "Ray/Intersection.hpp"
#include "Records.h"
#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "scene.hpp"
#include "Integrator/TraceHelper.h"
enum class VertexType {
    Camera,
    Light,
    Surface,
    Medium
};
class Camera;
class Light;
class Medium;

struct PathState {
    PathState(Sampler& sampler) : sampler(sampler) {
        bounce = 0;
    }

    Sampler& sampler;
    int      bounce;
    Medium*  medium = nullptr;
};

class PathVertex {
public:
    VertexType type;
    Spectrum   beta;
    Float      pdfBack, pdfFwd;
    union VertexRecord {
        LightRecord   lightRecord;
        CameraRecord  cameraRecord;
        SurfaceRecord surfaceRecord;
        MediumRecord  mediumRecord;
        VertexRecord() {}
        VertexRecord(const SurfaceRecord& surfaceRecord) {
        }
        ~VertexRecord() {}
    };
    union VertexSampler {
        const Light*  light;
        const Camera* camera;
        const BSDF*   bsdf;
        const Medium* medium;
    };
    VertexRecord  _record;
    VertexSampler _sampler;

    bool diarc = false;

public:
    PathVertex() {}
    //    PathVertex(PathVertex & other) = delete;
    PathVertex(const Light* light, Spectrum Le, Float pdf) {
        type           = VertexType::Light;
        _sampler.light = light;
        beta           = Le;
        pdfFwd         = pdf;
    }
    void initLight(const Light* light, Float lightPdf) {
        type                         = VertexType::Light;
        _sampler.light               = light;
        _record.lightRecord.lightPdf = lightPdf;
    }

    PathVertex(const Light* light, Float lightPdf) {
        initLight(light, lightPdf);
    }
    void initCamera(const Camera* camera, vec2 point) {
        type                       = VertexType::Camera;
        _sampler.camera            = camera;
        _record.cameraRecord.pixel = point;
    }

    PathVertex(const Camera* camera, vec2 point) {
        initCamera(camera, point);
    }

    PathVertex(const Camera* camera, PositionAndDirectionSample sample) {
        type                        = VertexType::Camera;
        _sampler.camera             = camera;
        _record.cameraRecord.sample = sample;
        beta                        = sample.weight / (sample.dirPdf * sample.posPdf);
    }

    PathVertex(const Light* light, PositionAndDirectionSample sample, Float lightPdf) {
        type           = VertexType::Light;
        _sampler.light = light;
        ;
        _record.lightRecord.sample   = sample;
        _record.lightRecord.lightPdf = lightPdf;
        beta                         = sample.weight / lightPdf;
        pdfFwd                       = 0;
    }

    //    void initSurface(const SurfaceRecord & record,const Spectrum  &beta){
    //        type = VertexType::Surface;
    //        _record = record;
    //        _sampler.bsdf = record.its.bsdf;
    //        this->beta = beta;
    //    }
    //As the tungsten author said,
    // this is an ugly situation,
    // and hopefully, it can be resolved during the refactoring of the surface event.
    void pointerFixUp() {
        _record.surfaceRecord.event.its = &_record.surfaceRecord.its;
    }
    PathVertex(const SurfaceRecord& record, const Spectrum& beta) : type(VertexType::Surface)

    {
        _record.surfaceRecord = record;
        _sampler.bsdf         = record.its.bsdf;
        this->beta            = beta;
    }

    bool sampleNext(const Scene& scene, bool adjoint, PathState& state, PathVertex* prev, PathVertex& next);
    bool sampleRootVertex(PathState& state);

    inline const Light* getLight() const {
        if (type == VertexType::Light)
            return _sampler.light;
        if (type == VertexType::Surface)
            return _record.surfaceRecord.its.primitive->getAreaLight();
        return nullptr;
    }

    inline bool isLight() const {
        return type == VertexType::Light || (isSurface() && _record.surfaceRecord.its.primitive->getAreaLight());
    }
    inline bool isCamera() const {
        return type == VertexType::Camera;
    }

    inline bool isSurface() const {
        return type == VertexType::Surface;
    }
    inline bool isMedium() const {
        return type == VertexType::Medium;
    }

    inline vec3 ng() const {
        if (isSurface())
            return _record.surfaceRecord.its.Ng;
        if (isCamera())
            return _record.cameraRecord.sample.n;
        if (isLight())
            return _record.lightRecord.sample.n;
    }

    Float cosFactor(const PathVertex& v) const {
        if (isMedium())
            return 1.f;
        auto dir    = normalize(v.pos() - pos());
        auto result = absDot(v.ng(), dir);
        if (result == 0) {
            DebugBreak();
        }
        return result;
    }

    inline bool isDelta() const {
        return diarc;
    }
    Float ConvertDensity(Float pdf, const PathVertex& next) const {
        // Return solid angle density if _next_ is an infinite area light
        if (next.isInfiniteLight()) return pdf;
        vec3 w = next.pos() - pos();
        if (length2(w) == 0) return 0;
        Float invDist2 = 1 / length2(w);
        if (next.isSurface())
            pdf *= absDot(next.ng(), w * std::sqrt(invDist2));
        return pdf * invDist2;
    }

    bool canConnect() const;

    vec3 pos() const;

    Spectrum eval(const PathVertex& vertex, bool adjoint) const;

    //对于相机和光源节点 prev可能为nullptr
    Float pdf(const PathVertex* prev, const PathVertex& next) const;
    //计算作为光源，给定下一个节点v的pdf
    Float pdfLight(const PathVertex& v) const;
    //计算被选中光源，采样出给定方向的pdf
    Float pdfLightOrigin(const PathVertex&     v,
                         const Distribution1D& lightDistr,
                         const std::map<const Light*, size_t>) const;

    bool isInfiniteLight() const {
        return false;
    }
};