//2022/9/4
#pragma  once

#include "spdlog/spdlog-inl.h"
#include <embree3/rtcore.h>
#include <optional>

#define   _CHECK_RTC_DEVICE  assert(rtcGetDeviceError(EmbreeUtils::getDevice())==RTC_ERROR_NONE);
#define   _LOG_RTC_DEVICE_STATES_IF_ERROR if(rtcGetDeviceError(EmbreeUtils::getDevice())!=RTC_ERROR_NONE) \
                                            spdlog::info(rtcGetDeviceError(EmbreeUtils::getDevice()));

class Intersection;
class Ray;
namespace  EmbreeUtils {
    struct RTCRayHit1 : RTCRayHit{
        Intersection * its = nullptr ;
    };

    inline RTCRay *  RTCRay_(RTCRayN *  ray) {
        return reinterpret_cast<RTCRay*>(ray); }

    inline RTCRayHit* RTCRayHit_(RTCRayHitN * rayHit) {
    return reinterpret_cast<RTCRayHit*>(rayHit); }

    inline RTCRayHit1* RTCRayHit1_(RTCRayHitN * rayHit) {
    return reinterpret_cast<RTCRayHit1*>(rayHit); }

    inline RTCRayHit1* RTCRayHit1_(RTCRayHit * rayHit) {
        return reinterpret_cast<RTCRayHit1*>(rayHit); }



    RTCDevice getDevice();
    void initDevice();

    void convertRay(const Ray *ray,RTCRay * rayhit );
    void convertRay(const Ray *ray,RTCRayHit * rayhit );
    void convertRay(const RTCRay * rtcray,Ray * ray);
    void instanceBoundsFunc(const struct RTCBoundsFunctionArguments* args);
    void instanceIntersectFunc(const RTCIntersectFunctionNArguments* args);
    void instanceOccludedFunc(const RTCOccludedFunctionNArguments* args);
}
