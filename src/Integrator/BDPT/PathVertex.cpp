#include "PathVertex.h"
#include "Lights/Light.hpp"
#include "Bsdfs/Reflection.hpp"
#include "Camera/Camera.hpp"

bool PathVertex::canConnect() const {
    switch (type) {
        case VertexType::Medium:
            return true;
        case VertexType::Light:
            // return (sampler.light->flags )
            return (_sampler.light->flags & static_cast<int>(DeltaDirection)) == 0;
        case VertexType::Camera:
            return true;
        case VertexType::Surface:
            return !_sampler.bsdf->Pure(BSDF_PURE_SPECULR);
    }
    //   LOG(FATAL) << "Unhandled vertex type in IsConnectable()";
    return false;// NOTREACHED
}

vec3 PathVertex::pos() const {
    switch (type) {
        case VertexType::Medium:
            //return true;
        case VertexType::Light:
            // return (sampler.light->flags )
            return _record.lightRecord.sample.ray.o;
            //  return (sampler.light->flags & (int) LightFlags::DeltaDirection) == 0;
        case VertexType::Camera:
            return _record.cameraRecord.sample.ray.o;
        case VertexType::Surface:
            return _record.surfaceRecord.its.p;
    }
}

Spectrum PathVertex::eval(const PathVertex& vertex, bool adjoint) const {
    vec3 d = normalize(vertex.pos() - pos());
    switch (type) {
        case VertexType::Surface: {
            const auto& event = _record.surfaceRecord.event;
            return _sampler.bsdf->f(event.makeWarpQuery(event.toLocal(d), event.wo), adjoint);
        }
        default:
            return Spectrum();
    }
    return Spectrum();
}

Float PathVertex::pdf(const PathVertex* prev, const PathVertex& next) const {

    //前一个节点是光源或者相机时需要特殊处理
    //pdf需要从角度空间转成面积空间
    if (isLight()) {
        //光源有两种处理 第一种是自己就是光源 第二种是前面一个是光源
        return pdfLight(next);
    }
    Float pdf = 0;

    vec3 wn = next.pos() - pos();
    vec3 wp;
    if (prev) {
        wp = normalize(prev->pos() - pos());
    } else {
        CHECK(isCamera(), "Not camera vertex but without prev vertex");
    }
    if (isSurface()) {
        wn  = normalize(wn);
        wn  = _record.surfaceRecord.event.toLocal(wn);
        wp  = _record.surfaceRecord.event.toLocal(wp);
        pdf = _sampler.bsdf->Pdf(_record.surfaceRecord.event.makeWarpQuery(wn, wp));
    }
    if (isCamera())
        _sampler.camera->pdfRay(Ray(_record.cameraRecord.sample.ray.o, normalize(wn)), nullptr, &pdf);
    if (isMedium())
        TODO("Impl medium pdf");
    //   pdf = 1;
    //return cosFactor(next);
    //return distance2(pos(),next.pos());
    return pdf * cosFactor(next) / distance2(pos(), next.pos());
}

Float PathVertex::pdfLight(const PathVertex& v) const {
    vec3  d        = v.pos() - pos();
    Float invDist2 = 1 / length2(d);
    d *= invDist2;
    Float pdf = 0;
    if (isInfiniteLight()) {
        TODO("Handle Infinite");
    } else {
        assert(isLight());
        const Light* light = getLight();
        // Compute sampling density for non-infinite light sources
        Float pdfDir;
        light->pdfDirect(Ray(pos(), d), ng(), nullptr, &pdfDir);
        pdf = pdfDir * invDist2;
        if (v.isSurface())
            pdf *= absDot(v.ng(), d);
    }
    return pdf;
}

Float PathVertex::pdfLightOrigin(const PathVertex& v, const Distribution1D& lightDistr, const std::map<const Light*, size_t> map) const {
    vec3  d        = v.pos() - pos();
    Float invDist2 = 1 / length2(d);
    d *= invDist2;
    Float pdf = 0;
    if (isInfiniteLight()) {
        TODO("Handle Infinite");
    } else {
        const Light* light = getLight();
        // Compute sampling density for non-infinite light sources
        Float pdfPos;
        light->pdfDirect(Ray(pos(), d), ng(), &pdfPos, nullptr);
        pdf = pdfPos * lightDistr.DiscretePDF(0);
    }
    return pdf;
}

bool PathVertex::sampleNext(const Scene& scene, bool adjoint, PathState& state, PathVertex* prev, PathVertex& next) {
    Spectrum weight(1);
    Float    pdf(1);
    Ray      ray;
    switch (type) {
        case VertexType::Light: {
            auto& record = _record.lightRecord;
            ray          = record.sample.ray;
            pdf          = record.sample.dirPdf;
            weight       = Spectrum(1.f) / (record.sample.dirPdf);
            break;
        }
        case VertexType::Camera: {
            auto& record = _record.cameraRecord;
            pdf          = record.sample.dirPdf;
            ray          = record.sample.ray;
            break;
        }
        case VertexType::Medium: {
            break;
        }
        case VertexType::Surface: {
            auto&         record = _record.surfaceRecord;
            SurfaceEvent& event  = record.event;
            auto&         bsdf   = _sampler.bsdf;
            weight               = bsdf->sampleF(event, state.sampler.getNext2D(), adjoint);
            if (isBlack(weight) || event.pdf == 0) {
                return false;
            }
            if (event.sampleType & BSDF_SPECULAR)
                diarc = true;
            weight /= event.pdf;
            if (russian(state.bounce, state.sampler, weight * beta))
                return false;
            pdf           = event.pdf;
            ray           = event.sctterRay();
            prev->pdfBack = bsdf->Pdf(event.makeFlipQuery());
            break;
        }
    }
    auto its = scene.intersect(ray);
    if (!its) {
        return false;
    }
    SurfaceRecord record;
    record.its   = its.value();
    auto event   = makeLocalScatterEvent(&record.its);
    record.event = event;

    next = PathVertex(record, beta * weight);
    next.pointerFixUp();
    next.pdfFwd = pdf;
    state.bounce++;

    return true;
}

///Maybe it's better to separate the sampling process of pos and dir here.
bool PathVertex::sampleRootVertex(PathState& state) {

    switch (type) {
        case (VertexType::Camera): {
            auto  camera  = _sampler.camera;
            auto& record  = _record.cameraRecord;
            record.sample = camera->sampleRay(record.pixel, state.sampler.getNext2D(), state.sampler.getNext2D());
            beta          = record.sample.weight;
            pdfFwd        = record.sample.posPdf;
            return true;
        }
        case (VertexType::Light): {
            auto  light   = _sampler.light;
            auto& record  = _record.lightRecord;
            record.sample = light->sampleDirect(state.sampler.getNext2D(), state.sampler.getNext2D());
            if (record.sample.dirPdf == 0 || record.sample.posPdf == 0)
                return false;
            if (isInfiniteLight()) {
                //todo
            } else {
                beta   = record.sample.weight * absDot(record.sample.n, record.sample.ray.d) / record.lightPdf / record.sample.posPdf;
                pdfFwd = record.sample.posPdf * record.lightPdf;
            }
            return true;
        }
        default:
            return false;
    }
    return false;
}





//PathVertex::VertexRecord::VertexRecord() {
//  
//}