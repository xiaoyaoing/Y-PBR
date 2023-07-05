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
            return (_sampler.light->flags & (int) LightFlags::DeltaDirection) == 0;
        case VertexType::Camera:
            return true;
        case VertexType::Surface:
            return !_sampler.bsdf->Pure(BSDF_PURE_SPECULR);
    }
    //   LOG(FATAL) << "Unhandled vertex type in IsConnectable()";
    return false;  // NOTREACHED

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

Spectrum PathVertex::eval(const PathVertex &vertex, bool adjoint) const {
    vec3 d = normalize(vertex.pos() - pos());
    switch (type) {
        case VertexType::Surface: {
            const auto &event = _record.surfaceRecord.event;
          //  return (1.f+event.toLocal(d))/2.f;
            return _sampler.bsdf->f(event.makeWarpQuery( event.wi,event.toLocal(d)), false);
        }
        default:
            return Spectrum();
    }
    return Spectrum();
}

bool PathVertex::sampleNext(const Scene &scene, bool adjoint, PathState &state, PathVertex *prev, PathVertex &next) {
    Spectrum weight(1);
    Float pdf(1);
    Ray ray;
    switch (type) {
        case VertexType::Light: {
            auto &record = _record.lightRecord;
          //  pdf = record.sample.dirPdf;
//           weight = record.sample.radiance;
            ray = record.sample.ray;
            pdf = record.sample.dirPdf ;
            weight = Spectrum (dot(record.sample.n,ray.d) /(record.sample.dirPdf)) ;
            break;
        }
        case VertexType::Camera: {
            auto &record = _record.cameraRecord;
//            weight = record.sample.radiance;
            pdf = record.sample.dirPdf;
            ray = record.sample.ray;
            break;
        }
        case VertexType::Medium: {
            break;
        }
        case VertexType::Surface: {
            auto &record = _record.surfaceRecord;
            SurfaceEvent &event = record.event;
            auto &bsdf = event.its->bsdf;
            weight = bsdf->sampleF(event, state.sampler.getNext2D(), adjoint);
            if(isBlack(weight) || event.pdf == 0){
                return false;
            }
            weight/=event.pdf;
            if(isinf(weight[0])){
                int k = 1;
            }
            if (russian(state.bounce, state.sampler,  weight * beta))
                return false;
            pdf = event.pdf;
            ray = event.sctterRay();
            prev->pdfBack = bsdf->Pdf(event.makeFlipQuery());
            break;
        }
    }
    auto its = scene.intersect(ray);
    if (!its) {
        return false;
    }
    SurfaceRecord record;
    record.its = its.value();
    auto event = makeLocalScatterEvent(&record.its);
    record.event = event;
    if(hasNan(weight) || isinf(weight[0])){
        int k = 1;
    }
    next._record = VertexRecord();
    next._record.surfaceRecord = SurfaceRecord();
    next = PathVertex(record, beta * weight);
    next.pdfFwd = pdf;
    state.bounce++;
    return true;
}


///Maybe it's better to separate the sampling process of pos and dir here.
bool PathVertex::sampleRootVertex(PathState &state) {

    switch (type) {
        case (VertexType::Camera): {
            auto camera = _sampler.camera;
            auto &record = _record.cameraRecord;
            record.sample = camera->sampleRay(record.pixel, state.sampler.getNext2D(), state.sampler.getNext2D());
            beta = record.sample.weight;
            pdfFwd = record.sample.posPdf;
            return true;
        }
        case (VertexType::Light): {
            auto light = _sampler.light;
            auto &record = _record.lightRecord;
            record.sample = light->sampleDirect(state.sampler.getNext2D(), state.sampler.getNext2D());
            if (isInfiniteLight()) {
                //todo
            } else {
                beta = record.sample.weight/ record.lightPdf / record.sample.posPdf;
                if(beta[0]==0){
                    int k = 1;
                }
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

