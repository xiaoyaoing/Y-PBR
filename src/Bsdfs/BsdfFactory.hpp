#include "Reflection.hpp"
#include "Conductor.hpp"

namespace BSDFFactory{

std::shared_ptr<BSDF> LoadBsdfFromJson(Json j);

std::unordered_map<std::string,std::shared_ptr<BSDF>>
                      LoadBsdfsFromJson(Json j);

}
