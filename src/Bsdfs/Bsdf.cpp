#include "Bsdf.hpp"


void from_json(const nlohmann::json &j, Bsdf &m){

}




void std::from_json(const nlohmann::json &j, std::shared_ptr<Bsdf> &m)
{
    m = std::make_shared<Bsdf>(j.get<Bsdf>());
}

Spectrum Bsdf::f(const vec3 & wo, const vec3 & wi) const {
    return Spectrum();
}
