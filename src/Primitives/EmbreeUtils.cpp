#include "EmbreeUtils.hpp"
#include "Primitive.hpp"
#include "Ray/Intersection.hpp"
#include "spdlog/spdlog.h"

namespace EmbreeUtils {

    static RTCDevice globalDevice = nullptr;

    struct RTCIntersectContext1 : RTCIntersectContext {
        Intersection * its;
    };


    RTCDevice getDevice( ) {
        if ( globalDevice == nullptr )
            initDevice();
        return globalDevice;
    }

    void initDevice( ) {
        globalDevice = rtcNewDevice(nullptr);
        RTCError error = rtcGetDeviceError(globalDevice);
        if ( error != RTC_ERROR_NONE ) {
            std::string err_str("Embree device creation error: ");
            err_str.append(std::to_string(error));
            throw std::runtime_error(err_str);
        }
    }

    void convertRay(const RTCRay * rtcray, Ray * ray) {
        ray->o = vec3(rtcray->org_x, rtcray->org_y, rtcray->org_z);
        ray->d = vec3(rtcray->dir_x, rtcray->dir_y, rtcray->dir_z);
        ray->nearT = rtcray->tnear;
        ray->farT = rtcray->tfar;
    }

    void convertRay(const Ray * ray, RTCRay * rtcRay) {
        rtcRay->org_x = ray->o.x;
        rtcRay->org_y = ray->o.y;
        rtcRay->org_z = ray->o.z;
        rtcRay->dir_x = ray->d.x;
        rtcRay->dir_y = ray->d.y;
        rtcRay->dir_z = ray->d.z;
        rtcRay->tnear = ray->nearT;
        rtcRay->tfar = ray->farT;
    }


    void convertRay(const Ray * ray, RTCRayHit * rayhit) {
        rayhit->ray.org_x = ray->o.x;
        rayhit->ray.org_y = ray->o.y;
        rayhit->ray.org_z = ray->o.z;
        rayhit->ray.dir_x = ray->d.x;
        rayhit->ray.dir_y = ray->d.y;
        rayhit->ray.dir_z = ray->d.z;
        rayhit->ray.tnear = ray->nearT;
        rayhit->ray.tfar = ray->farT;
        rayhit->hit.geomID = RTC_INVALID_GEOMETRY_ID;
    }

    void instanceBoundsFunc(const struct RTCBoundsFunctionArguments * args) {
        const Primitive * instance = (const Primitive *) args->geometryUserPtr;
        RTCBounds * bounds_o = args->bounds_o;
        vec3 bMin = instance->BB().pMin;
        vec3 bMax = instance->BB().pMax;
        bounds_o->lower_x = bMin.x;
        bounds_o->lower_y = bMin.y;
        bounds_o->lower_z = bMin.z;
        bounds_o->align0 = 0;
        bounds_o->upper_x = bMax.x;
        bounds_o->upper_y = bMax.y;
        bounds_o->upper_z = bMax.z;
        bounds_o->align1 = 0;
    }

    void instanceIntersectFunc(const RTCIntersectFunctionNArguments * args) {
        //spdlog::info("intersectFunc");
        void * ptr = args->geometryUserPtr;
        Primitive * primitive = (Primitive *) ptr;
        RTCRayHit1 * hit = RTCRayHit1_(args->rayhit);
        Ray ray;
        convertRay(& hit->ray, & ray);
        std::optional < Intersection > its = primitive->intersect(ray);
        if ( its.has_value() ) {
            hit->hit.geomID = args->geomID;
            hit->its = & ( its.value() );
            hit->ray.tfar = ray.farT;
        }
    }


    void instanceOccludedFunc(const RTCOccludedFunctionNArguments * args) {
        void * ptr = args->geometryUserPtr;
        Primitive * primitive = (Primitive *) ptr;
        RTCRay * rtcRay = RTCRay_(args->ray);
        Ray ray;
        convertRay(rtcRay, & ray);
        if (primitive->occluded(ray) ) {
            RTCRay_(args->ray)->tfar = _NEG_INFINY;
        }
    }
}