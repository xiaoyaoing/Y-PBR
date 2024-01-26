#include "Common/math.hpp"

namespace Fresnel {

    static inline Float dielectricReflectance(Float eta, Float cosThetaI, Float& cosThetaT) {
        if (cosThetaI < 0.0) {
            eta       = 1.0 / eta;
            cosThetaI = -cosThetaI;
        }
        Float sinThetaTSq = eta * eta * (1.0f - cosThetaI * cosThetaI);
        if (sinThetaTSq > 1.0) {
            cosThetaT = 0.0;
            return 1.0;
        }
        cosThetaT = std::sqrt(std::max(1.0 - sinThetaTSq, 0.0));

        Float Rs = (eta * cosThetaI - cosThetaT) / (eta * cosThetaI + cosThetaT);
        Float Rp = (eta * cosThetaT - cosThetaI) / (eta * cosThetaT + cosThetaI);

        return (Rs * Rs + Rp * Rp) * 0.5;
    }

    static inline Float dielectricReflectance(Float eta, Float cosThetaI) {
        Float cosThetaT;
        return dielectricReflectance(eta, cosThetaI, cosThetaT);
    }

    static inline Float SchlickApproxFresnel(Float eta, Float cosTheta) {
        float r0 = (eta * eta - 2 * eta + 1) / (eta * eta + 2 * eta + 1);
        return r0 + (1 - r0) * pow(1 - cosTheta, 5);
    }

    static inline Float conductorReflectance(Float eta, Float k, Float cosThetaI) {
        Float cosThetaISq = cosThetaI * cosThetaI;
        Float sinThetaISq = std::max(1.0f - cosThetaISq, 0.0f);
        Float sinThetaIQu = sinThetaISq * sinThetaISq;

        Float innerTerm  = eta * eta - k * k - sinThetaISq;
        Float aSqPlusBSq = std::sqrt(std::max(innerTerm * innerTerm + 4.0f * eta * eta * k * k, 0.0f));
        Float a          = std::sqrt(std::max((aSqPlusBSq + innerTerm) * 0.5f, 0.0f));

        Float Rs = ((aSqPlusBSq + cosThetaISq) - (2.0f * a * cosThetaI)) /
                   ((aSqPlusBSq + cosThetaISq) + (2.0f * a * cosThetaI));
        Float Rp = ((cosThetaISq * aSqPlusBSq + sinThetaIQu) - (2.0f * a * cosThetaI * sinThetaISq)) /
                   ((cosThetaISq * aSqPlusBSq + sinThetaIQu) + (2.0f * a * cosThetaI * sinThetaISq));

        return 0.5f * (Rs + Rs * Rp);
    }

    static inline vec3 conductorReflectance(const vec3& eta, const vec3& k, float cosThetaI) {
        return vec3(
            conductorReflectance(eta.x, k.x, cosThetaI),
            conductorReflectance(eta.y, k.y, cosThetaI),
            conductorReflectance(eta.z, k.z, cosThetaI));
    }

    static inline Float diffuseReflectance(Float eta, Float sampleCount) {
        double diffuseFresnel = 0.0;
        float  fb             = Fresnel::dielectricReflectance(eta, 0.0f);
        for (int i = 1; i <= sampleCount; ++i) {
            float cosThetaSq = float(i) / sampleCount;
            float fa         = Fresnel::dielectricReflectance(eta, std::min(std::sqrt(cosThetaSq), 1.0f));
            diffuseFresnel += double(fa + fb) * (0.5 / sampleCount);
            fb = fa;
        }

        return float(diffuseFresnel);
    }
}