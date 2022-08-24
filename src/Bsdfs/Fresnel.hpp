#include "Common/math.hpp"

namespace Fresnel {

static inline  Float DielectricReflectance(Float eta,Float cosThetaI,Float cosThetaT){
    if (cosThetaI < 0.0) {
        eta = 1.0/eta;
        cosThetaI = -cosThetaI;
    }
    float sinThetaTSq = eta*eta*(1.0f - cosThetaI*cosThetaI);
    if (sinThetaTSq > 1.0) {
        cosThetaT = 0.0;
        return 1.0;
    }
    cosThetaT = std::sqrt(std::max(1.0 - sinThetaTSq, 0.0));

    float Rs = (eta*cosThetaI - cosThetaT)/(eta*cosThetaI + cosThetaT);
    float Rp = (eta*cosThetaT - cosThetaI)/(eta*cosThetaT + cosThetaI);

    return (Rs*Rs + Rp*Rp)*0.5;
}

}
