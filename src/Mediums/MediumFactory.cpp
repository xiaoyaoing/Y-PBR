#include "MediumFactory.hpp"
#include "Medium.hpp"

template < class T >
std::shared_ptr < T > LoadSimpleHelper(const Json & json) { return std::make_shared < T >(json); }

std::unordered_map < std::string, std::function < std::shared_ptr < Medium >(const Json &)>> loadMediumMap = {
        {"homogeneous", LoadSimpleHelper < Homogeneous >},
};

std::unordered_map < std::string, std::shared_ptr < Medium>> MediumFactory::loadMediumsFromJson(const Json & json) {
    std::unordered_map<std::string, std::shared_ptr < Medium>> mediums;
    for ( const Json & mediumJson: json ) {
        std::string mediumType = getOptional(mediumJson, "type", std::string("homogenous"));
        std::string mediumName = mediumJson["name"];
        mediums[mediumName] = (loadMediumMap[mediumType](mediumJson));
    }
    return mediums;
}
