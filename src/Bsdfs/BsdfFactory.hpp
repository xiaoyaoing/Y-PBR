#include "Reflection.hpp"

namespace BsdfFactory{

std::shared_ptr<Bsdf> LoadBsdfFromJson(nlohmann::json j);

std::unordered_map<std::string,std::shared_ptr<Bsdf>>
                      LoadBsdfsFromJson(nlohmann::json j);

}
