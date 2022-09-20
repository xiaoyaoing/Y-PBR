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

    enum SplitMethod { SAH, HLBVH, Middle, EqualCounts,PBRTSAH };

public:
    BVHAccel( //const BoundingBox &BB,
        std::vector < std::shared_ptr < Primitive>>   p,
        const SplitMethod splitMethod = PBRTSAH,
        const uint32 maxPrimInNode = 4
        );

    std::optional<Intersection> intersect(const Ray & ray) const ;

    bool intersectP(const Ray & ray) const ;

private:
    int FlattenTree(BuildNode *node, int *offset);

    std::vector<std::shared_ptr<Primitive>> primitives;
    const int maxPrimNode;
    SplitMethod splitMethod;
    LinearNode *nodes = nullptr;
    BuildNode * root;
};

std::unique_ptr<BVHAccel> CreateBVH(const nlohmann::json &j,
                                    const std::vector<std::shared_ptr<Primitive>> &
                                    );