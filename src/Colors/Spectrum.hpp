#include "../Common/math.hpp"
#pragma  once
using Spectrum = vec3;

inline bool isBlack(const Spectrum & color ){
    return color.x==0 && color.y==0 || color.z==0;
}

//void from_json(const nlohmann::json &j, Spectrum & spectrum);
