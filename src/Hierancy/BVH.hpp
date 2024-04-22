#pragma once
#include "../Primitives/Primitive.hpp"

struct alignas(32) LinearNode {
    Bounds3 bounds;
    union {
        int primitivesOffset; // leaf
        int secondChildOffset;// interior
    };
    uint16_t nPrimitives;// 0 -> interior node
    uint8_t  axis;
};

struct BuildNode;
struct BVHPrimitiveInfo;
//BVHAccel Mainly  Borrowed from PBRTv3
//Not really Used in Y_PBR
//Embree is much faster
class BVHAccel {

    BuildNode* RecursiveBuild(std::vector<BVHPrimitiveInfo>& primitiveInfo, int start, int end, std::vector<std::shared_ptr<Primitive>>& orderedPrims);

public:
    enum SplitMethod { SAH,
                       HLBVH,
                       Middle,
                       EqualCounts,
                       PBRTSAH };
    BVHAccel(//const BoundingBox &BB,
        std::vector<std::shared_ptr<Primitive>> p,
        const SplitMethod                       splitMethod   = PBRTSAH,
        const uint32                            maxPrimInNode = 4);

    std::optional<Intersection> intersect(const Ray& ray) const;
    bool                        intersectP(const Ray& ray) const;

    template<typename LAMBDA>
    void trace(Ray& ray, LAMBDA intersector) const {
        if (!nodes) {
            return;
        }
        //  auto temp=root->intersect(primitives,_ray);
        //  LOGI("Search Count {0}",IntersectionCount);
        //IntersectionCount = 0;
        //  return temp;

        vec3 invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
        int  dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};

        int                         toVisit[64];
        int                         visitIdx   = 0;
        int                         curNodeIdx = 0;
        std::optional<Intersection> res;
        int                         searchCount = 0;
        while (true) {
            //LOGI("CurNode {0}",curNodeIdx);
            LinearNode& curNode = nodes[curNodeIdx];
            searchCount++;
            if (curNode.bounds.IntersectP(ray, invDir, dirIsNeg)) {
                // node case
                if (curNode.nPrimitives > 0) {
                    //LOGI("LeafNode {0} ",curNodeIdx);
                    for (int i = 0; i < curNode.nPrimitives; ++i) {
                        intersector(ray, primitives[curNode.primitivesOffset + i]->primId);
                    }
                    if (visitIdx == 0) break;//no node to visit
                    curNodeIdx = toVisit[--visitIdx];
                }
                //interior node case
                else {
                    //near
                    if (dirIsNeg[curNode.axis]) {
                        //LOGI("push {0} to visitList",curNode.secondChildOffset);
                        curNodeIdx += 1;
                        toVisit[visitIdx++] = curNode.secondChildOffset;
                    } else {
                        toVisit[visitIdx++] = curNodeIdx + 1;
                        curNodeIdx          = curNode.secondChildOffset;
                        //LOGI("push {0} to visitList",curNodeIdx+1);
                    }
                }
            } else {
                //LOGI("CurNode not hit  {0}",curNodeIdx);
                if (visitIdx == 0) break;
                curNodeIdx = toVisit[--visitIdx];
            }
        }
    }

private:
    int FlattenTree(BuildNode* node, int* offset);

    std::vector<std::shared_ptr<Primitive>> primitives;
    const int                               maxPrimNode;
    SplitMethod                             splitMethod;
    LinearNode*                             nodes = nullptr;
    BuildNode*                              root;
};

std::unique_ptr<BVHAccel> CreateBVH(const Json& j,
                                    const std::vector<std::shared_ptr<Primitive>>&);