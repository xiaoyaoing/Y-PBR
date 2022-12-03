#include "PhaseFunctionFactory.hpp"
#include "PhaseFunction.hpp"
namespace PhaseFunctionFactory {

    template < class T >
    std::shared_ptr < T > LoadSimpleHelper(const Json & json) { return std::make_shared < T >(json); }

    std::unordered_map < std::string, std::function < std::shared_ptr < PhaseFunction >(const Json &)>> loadPhaseFunctionMap = {
            {"isotropic", LoadSimpleHelper <IsotropicPhaseFunction >},
    };
    std::shared_ptr < PhaseFunction > loadPhaseFromJson(const Json & json){
        std::string phaseType = json["type"];
        return loadPhaseFunctionMap[phaseType](json);
    }
}