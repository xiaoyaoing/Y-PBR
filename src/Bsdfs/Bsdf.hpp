#pragma  once

#include "nlohmann/json.hpp"

class Bsdf{

};

void from_json(const nlohmann::json &j, Bsdf &m);

namespace std
{
    void from_json(const nlohmann::json &j, shared_ptr<Bsdf> &m);
}