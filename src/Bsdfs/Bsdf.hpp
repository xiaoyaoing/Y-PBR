#pragma  once

#include "nlohmann/json.hpp"
#include "../Colors/Spectrum.hpp"

class Bsdf{
public:
    virtual Spectrum f(const vec3 & wo,const vec3 & wi) const ;
};

void from_json(const nlohmann::json &j, Bsdf &m);

namespace std
{
    void from_json(const nlohmann::json &j, shared_ptr<Bsdf> &m);
}