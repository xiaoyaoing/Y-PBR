#include "Common/Json.hpp"
class Medium;
namespace MediumFactory{
    std::unordered_map < std::string, std::shared_ptr < Medium>> loadMediumsFromJson(const Json & json);
}