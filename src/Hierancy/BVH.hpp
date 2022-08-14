#pragma  once
#include "../Primitives/Primitive.hpp"

struct LinearNode;
struct BuildNode;
struct BVHPrimitiveInfo;


//BVHAccel Mainly  Borrowed from PBRTv3


class BVHAccel
{

    BuildNode * RecursiveBuild( std::vector<BVHPrimitiveInfo> &primitiveInfo, int start, int end,
                                std::vector<std::shared_ptr<Primitive>> &orderedPrims );

    std::optional<Intersection> intersect(const Ray & ray);

    bool intersectP(const Ray & ray);

    enum SplitMethod { SAH, HLBVH, Middle, EqualCounts };

public:
    BVHAccel( //const BoundingBox &BB,
        const std::vector<std::shared_ptr<Primitive>> & primitives,
        const SplitMethod splitMethod = Middle,
        const uint32 maxPrimInNode = 4
        );

private:
    int FlattenTree(BuildNode *node, int *offset);

    std::vector<std::shared_ptr<Primitive>> primitives;
    SplitMethod splitMethod;
    std::vector<LinearNode> nodes;
};

std::shared_ptr<BVHAccel> CreateBVH(const nlohmann::json &j,
                                    const std::vector<std::shared_ptr<Primitive>>
                                    );