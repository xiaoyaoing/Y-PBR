#include "Reflection.hpp"
#include "Conductor.hpp"
#include "BSSRDF.hpp"
namespace BSDFFactory {

    std::shared_ptr<BSDF> LoadBsdfFromJson(const Json& j);

    std::unordered_map<std::string, std::shared_ptr<BSDF>>
    LoadBsdfsFromJson(const Json& j);

    std::unordered_map<std::string, std::shared_ptr<BSSRDF>>
    LoadBssrdfsFromJson(const Json& j);
}