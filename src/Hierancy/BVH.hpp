#pragma  once
#include "../Primitives/Primitive.hpp"


class BVHAccel
{

    struct BuildNode
    {
        BuildNode() { }

        bool leaf()
        {
            return children.empty();
        }

        BoundingBox BB;
        std::vector<std::shared_ptr<BuildNode>> children;
        std::vector<std::shared_ptr<Primitive>> primitives;
        uint32 df_idx; // depth-first index in tree
    };

    enum SplitMethod { SAH, HLBVH, Middle, EqualCounts };

public:
    BVHAccel(const BoundingBox &BB,
        const std::vector<std::shared_ptr<Primitive>> & primitives,
        const SplitMethod splitMethod = Middle,
        const uint32 maxPrimInNode = 4
        );
};

std::shared_ptr<BVHAccel> CreateBVH(const nlohmann::json &j,
                                    const std::vector<std::shared_ptr<Primitive>>
                                    );