#include "Common/math.hpp"
#include <vector>

namespace CurveIO {
    struct CurveData {
        std::vector<uint32>* curveEnds  = nullptr;
        std::vector<vec4>*   nodeData   = nullptr;
        std::vector<vec3>*   nodeColor  = nullptr;
        std::vector<vec3>*   nodeNormal = nullptr;
    };

    bool load(const std::string& path, CurveData& data);
    bool save(const std::string& path, const CurveData& data);

}