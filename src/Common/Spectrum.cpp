#include "Spectrum.hpp"

bool hasNan(const Spectrum& color) {

    bool hasNan = isnan(color[0]) || isnan(color[1]) || isnan(color[2]);
    if (hasNan)
        throw std::runtime_error("Nan Error");
    return hasNan;
}