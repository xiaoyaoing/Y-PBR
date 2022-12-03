#include "Common/Json.hpp"

class PhaseFunction;
namespace PhaseFunctionFactory {
    std::shared_ptr < PhaseFunction > loadPhaseFromJson(const Json & json);
}

