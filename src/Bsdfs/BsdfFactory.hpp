#include "Reflection.hpp"
#include "Conductor.hpp"

namespace BSDFFactory{

std::shared_ptr<BSDF> LoadBsdfFromJson(nlohmann::json j);

std::unordered_map<std::string,std::shared_ptr<BSDF>>
                      LoadBsdfsFromJson(nlohmann::json j);

}
