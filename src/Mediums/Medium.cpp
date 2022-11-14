#include "Medium.hpp"

Float Homogeneous::TR(const Ray & ray) const {
    return 1 * exp((ray.nearT - ray.farT) * sigmaT);
}
