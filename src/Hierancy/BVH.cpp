#include "BVH.hpp"

BVHAccel::BVHAccel(const std::vector < std::shared_ptr < Primitive>> & primitives, const SplitMethod splitMethod,
                   const uint32 maxPrimInNode) {
    if(primitives.empty()) return ;
}

std::shared_ptr<BVHAccel> CreateBVH(const nlohmann::json &j,
                                    const std::vector<std::shared_ptr<Primitive>> primitives
){

    return std::make_shared <BVHAccel>(primitives);
}