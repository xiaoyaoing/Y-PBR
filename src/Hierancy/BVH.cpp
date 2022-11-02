#include "BVH.hpp"
#include "../Common/Macro.hpp"
#include "set"
#include "iostream"
#include "spdlog/spdlog.h"
static int BuildNodeCount = 0;
static int IntersectionCount = 0;
struct BuildNode {
    BuildNode( ) {}

    std::optional<Intersection> intersect(  const   std::vector < std::shared_ptr < Primitive>> & primitives
    ,Ray & ray){
        IntersectionCount ++;
        std::optional<Intersection> val;
        if(nPrimitives){
            for(int i=firstPrimOffset;i<firstPrimOffset+nPrimitives;i++){
                auto its=primitives[i]->intersect(ray);
                if(its.has_value()) val=its;
            }
        }
        else {
        vec3 invDir(1/ray.d.x,1/ray.d.y,1/ray.d.z);
        int dirIsNeg[3] = {ray.d.x<0 ,ray.d.y<0 ,ray.d.z<0};
        if(!BB.IntersectP(ray,invDir,dirIsNeg))
            return std::nullopt;

        auto val1=children[0]->intersect(primitives,ray);
        if(val1.has_value())
            val=val1;
            val1=children[1]->intersect(primitives,ray);
        if(val1.has_value())
            val=val1;

        return val;
    }
    }

    void initLeaf(int off, int n, Bounds3 bounds3) {
        firstPrimOffset = off;
        nPrimitives = n;
        bounds3 = std::move(bounds3);
    }

    void initInterior(int dim, BuildNode * left, BuildNode * right) {
        nPrimitives = 0;
        splitAxis = dim;
        children[0] = left;
        children[1] = right;
    }

    //for debug
    void logInfo(bool logInterior,bool logLeaf){
        bool isLeaf=(nPrimitives>0);
        if(isLeaf && logLeaf){
            spdlog::info("Leaf Node start idx:{0} primNum:{1}",firstPrimOffset,nPrimitives);
        }
        if(!isLeaf && logInterior){
            spdlog::info("Interior node children : ");
        }

        if(!isLeaf){
            children[0]->logInfo(logInterior,logLeaf);
            children[1]->logInfo(logInterior,logLeaf);
        }
    }

    Bounds3 BB;
    BuildNode * children[2];
    int firstPrimOffset;
    int nPrimitives;
    int splitAxis;
};

struct BVHPrimitiveInfo {
    BVHPrimitiveInfo( ) {}

    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3 & bounds)
            : primitiveNumber(primitiveNumber),
              bounds(bounds),
              centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}

    static Float computeSurfaceArea(std::vector<BVHPrimitiveInfo> & infos,int start,int end){
        if(start==end) return 0;
        Bounds3 bounds3;
        for(int idx=start;idx<end;idx++){
            bounds3 = Union(bounds3,infos[idx].bounds);
        }
        return bounds3.SurfaceArea();
    }
    size_t primitiveNumber;
    Bounds3 bounds;
    vec3 centroid;
};


struct SAHBucketInfo {
    int count = 0;
    Bounds3 bounds;
};

struct alignas(32) LinearNode {
    Bounds3 bounds;
    union {
        int primitivesOffset;   // leaf
        int secondChildOffset;  // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;
};


BVHAccel::BVHAccel( std::vector < std::shared_ptr < Primitive>>   p,
                   const SplitMethod splitMethod,
                   const uint32 maxPrimInNode) : splitMethod(splitMethod),primitives(std::move(p)),maxPrimNode(maxPrimInNode)
                   {
    if ( primitives.empty() ) return;

    //Take the primitive and its spatial information apart8
    std::vector < BVHPrimitiveInfo > primitiveInfo(primitives.size());
    for ( size_t i = 0 ; i < primitives.size() ; ++ i )
        primitiveInfo[i] = {i, primitives[i]->BB()};

    std::vector < std::shared_ptr < Primitive>> orderedPrims;
    BuildNode * node = RecursiveBuild(primitiveInfo, 0, primitives.size(), orderedPrims);
    root = node;
    node->logInfo(true,true);
    spdlog::info("Build Num {0}", BuildNodeCount);
    primitives.swap(orderedPrims);
    nodes = new LinearNode[BuildNodeCount];
    int offset = 0;
    FlattenTree(node, & offset);
    spdlog::info("BVH Constructed");

    root=node;

#ifdef _DEBUG
    int leafNum = 0, interiorNum = 0, primitiveNum = 0;
    for ( int idx = 0 ; idx < BuildNodeCount ; idx ++ ) {
        LinearNode & node = nodes[idx];
        if ( node.nPrimitives ) {
            leafNum ++;
            primitiveNum += node.nPrimitives;
        } else {
            interiorNum ++;
        }
    }
    spdlog::info("NodeNum:{0} LeafNum:{1} interiorNum{2} PrimNum:{3}",
                 BuildNodeCount, leafNum, interiorNum, primitiveNum);
#endif
}

BuildNode * BVHAccel::RecursiveBuild(std::vector < BVHPrimitiveInfo > & primitiveInfo, int start, int end,
                                     std::vector < std::shared_ptr < Primitive>> & orderedPrims) {
    if ( start == end ) {
        return nullptr;
    }
    BuildNodeCount ++;
    Bounds3 bounds;
    for ( int i = start ; i < end ; i ++ ) {
        bounds = Union(bounds, primitiveInfo[i].bounds);
    }

    BuildNode * node = new BuildNode;
    node->BB = bounds;

    int nPrimitives = end - start;
    //Only single primitive
    if ( nPrimitives <= 10 ) {
        node->initLeaf(orderedPrims.size(), nPrimitives, bounds);
        for ( int i = start ; i < end ; ++ i ) {
            orderedPrims.push_back(primitives[primitiveInfo[i].primitiveNumber]);
        }
        return node;
    }

    Bounds3 centroidBounds;
    for ( int i = start ; i < end ; ++ i )
        centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
    int dim = centroidBounds.MaximumExtent();


    //All primitives are in one dimension
    if ( centroidBounds.pMax[dim] == centroidBounds.pMin[dim] ) {
        int primIdx = primitiveInfo[start].primitiveNumber;
        node->initLeaf(orderedPrims.size(), nPrimitives, bounds);
        for ( int i = start ; i < end ; ++ i ) {
            orderedPrims.push_back(primitives[primitiveInfo[i].primitiveNumber]);
        }
        return node;
    }

    int mid = ( start + end ) / 2;

    switch ( splitMethod ) {
        case Middle: {
            Float pmid = ( centroidBounds.pMin[dim] + centroidBounds.pMax[dim] ) / 2;
            BVHPrimitiveInfo * midPtr = std::partition(
                    & primitiveInfo[start], & primitiveInfo[end - 1] + 1,
                    [dim, pmid](const BVHPrimitiveInfo & pi) {
                        return pi.centroid[dim] < pmid;
                    });
            mid = midPtr - & primitiveInfo[0];
            if ( mid != start && mid != end ) break;
        }
        case EqualCounts: {
            mid = ( start + end ) / 2;
            std::nth_element(& primitiveInfo[start], & primitiveInfo[mid],
                             & primitiveInfo[end - 1] + 1,
                             [dim](const BVHPrimitiveInfo & a,
                                   const BVHPrimitiveInfo & b) {
                                 return a.centroid[dim] < b.centroid[dim];
                             });
        }
        case SAH : {
              int axis = centroidBounds.MaximumExtent();
              Float pMin = centroidBounds[0][dim];
              Float pMax = centroidBounds[1][dim];

              std::sort(primitiveInfo.begin()+start,primitiveInfo.begin()+end,
                        [dim](const BVHPrimitiveInfo & p1,const BVHPrimitiveInfo & p2){
                  return p1.centroid[dim] < p2.centroid[dim];
              });
              const int bucketNum = 20;

              Float minCost=std::numeric_limits<Float>::max();
              int bestSplitPos;
              for(int curBucket = 0; curBucket<bucketNum;curBucket++){
                  Float curP = pMin + curBucket * (pMax-pMin) / bucketNum;
                  BVHPrimitiveInfo tempInfo;
                  tempInfo.centroid[dim]  = curP;
                  int  splitPos = std::lower_bound(primitiveInfo.begin()+start,primitiveInfo.begin()+end,tempInfo,
                                             [&dim](const BVHPrimitiveInfo & info,const BVHPrimitiveInfo & tempInfo){
                      return info.centroid[dim]<tempInfo.centroid[dim];
                  }) - primitiveInfo.begin();

                  Float splitCost= BVHPrimitiveInfo::computeSurfaceArea(primitiveInfo,start,splitPos) * (splitPos-start)
                          + BVHPrimitiveInfo::computeSurfaceArea(primitiveInfo,splitPos,end) * (end-splitPos);
                  if(splitCost<minCost){
                      bestSplitPos = splitPos;
                      minCost = splitCost;
                  }
              }
              mid = bestSplitPos;
        }
            break;
        case PBRTSAH : {
            const int  nBuckets = 12;
            SAHBucketInfo buckets[nBuckets];
            for (int i = start; i < end; ++i) {
                int b = nBuckets *
                        centroidBounds.Offset(
                                primitiveInfo[i].centroid)[dim];
                if (b == nBuckets) b = nBuckets - 1;
                buckets[b].count++;
                buckets[b].bounds =
                        Union(buckets[b].bounds, primitiveInfo[i].bounds);
            }

            Float cost[nBuckets - 1];
            for (int i = 0; i < nBuckets - 1; ++i) {
                Bounds3 b0, b1;
                int count0 = 0, count1 = 0;
                for (int j = 0; j <= i; ++j) {
                    b0 = Union(b0, buckets[j].bounds);
                    count0 += buckets[j].count;
                }
                for (int j = i + 1; j < nBuckets; ++j) {
                    b1 = Union(b1, buckets[j].bounds);
                    count1 += buckets[j].count;
                }
                cost[i] = 1 +
                          (count0 * b0.SurfaceArea() +
                           count1 * b1.SurfaceArea()) /
                          bounds.SurfaceArea();
            }

            Float minCost = cost[0];
            int minCostSplitBucket = 0;
            for (int i = 1; i < nBuckets - 1; ++i) {
                if (cost[i] < minCost) {
                    minCost = cost[i];
                    minCostSplitBucket = i;
                }
            }

            Float leafCost = nPrimitives;
            if (nPrimitives > maxPrimNode || minCost < leafCost) {
                BVHPrimitiveInfo *pmid = std::partition(
                        &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                        [=](const BVHPrimitiveInfo &pi) {
                            int b = nBuckets *
                                    centroidBounds.Offset(pi.centroid)[dim];
                            if (b == nBuckets) b = nBuckets - 1;
                            return b <= minCostSplitBucket;
                        });
                mid = pmid - &primitiveInfo[0];
            } else {
                // Create leaf _BVHBuildNode_
                int firstPrimOffset = orderedPrims.size();
                for (int i = start; i < end; ++i) {
                    int primNum = primitiveInfo[i].primitiveNumber;
                    orderedPrims.push_back(primitives[primNum]);
                }
                node->initLeaf(firstPrimOffset, nPrimitives, bounds);
                return node;
            }

        }
        default:; //todo
    }

    node->initInterior(dim,
                       RecursiveBuild(primitiveInfo, start, mid, orderedPrims),
                       RecursiveBuild(primitiveInfo, mid, end, orderedPrims)
    );

    return node;
}

//To be honest, I don't know why PBRT flatten the tree, it's probably for performance reasons.
// Anyway, I'll follow this.
int BVHAccel::FlattenTree(BuildNode * node, int * offset) {
    LinearNode * linearNode = & nodes[* offset];
    linearNode->bounds = node->BB;
    int myOffset = ( * offset ) ++;
    if ( node->nPrimitives > 0 ) {
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        linearNode->nPrimitives = 0;
        linearNode->axis = node->splitAxis;
        FlattenTree(node->children[0], offset);
        linearNode->secondChildOffset = FlattenTree(node->children[1], offset);
    }
    return myOffset;
}

std::optional < Intersection > BVHAccel::intersect(const Ray & ray) const {
    if ( ! nodes ) {
        return std::nullopt;
    }
    Ray _ray(ray);
  //  auto temp=root->intersect(primitives,_ray);
  //  spdlog::info("Search Count {0}",IntersectionCount);
    //IntersectionCount = 0;
  //  return temp;


    vec3 invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};

    int toVisit[64];
    int visitIdx = 0;
    int curNodeIdx = 0;
    std::optional < Intersection > res;
    int searchCount = 0;
    while ( true ) {
        //spdlog::info("CurNode {0}",curNodeIdx);
        LinearNode & curNode = nodes[curNodeIdx];
        searchCount ++;
        if ( curNode.bounds.IntersectP(ray, invDir, dirIsNeg)   ) {
            // node case
            if ( curNode.nPrimitives > 0 ) {
                //spdlog::info("LeafNode {0} ",curNodeIdx);
                for ( int i = 0 ; i < curNode.nPrimitives ; ++ i ) {
                    auto its = primitives[curNode.primitivesOffset + i]->intersect(_ray);
                    if ( its.has_value() )
                        res = its;
                }
                if ( visitIdx == 0 ) break; //no node to visit
                curNodeIdx = toVisit[-- visitIdx];
            }
                //interior node case
            else {
                //near
                if ( dirIsNeg[curNode.axis] ) {
                    //spdlog::info("push {0} to visitList",curNode.secondChildOffset);
                    curNodeIdx += 1;
                    toVisit[visitIdx ++] = curNode.secondChildOffset;
                } else {
                    toVisit[visitIdx ++] = curNodeIdx + 1;
                    curNodeIdx = curNode.secondChildOffset;
                    //spdlog::info("push {0} to visitList",curNodeIdx+1);
                }
            }
        } else {
            //spdlog::info("CurNode not hit  {0}",curNodeIdx);
            if ( visitIdx == 0 ) break;
            curNodeIdx = toVisit[-- visitIdx];
        }
    }

     spdlog::info("SearchCount :{0}",searchCount);
    return res;

}

bool BVHAccel::intersectP(const Ray & ray) const {
    if ( ! nodes ) {
        return false;
    }
    Ray _ray(ray);
    //return root->intersect(primitives,_ray);

    vec3 invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};

    int toVisit[64];
    int visitIdx = 0;
    int curNodeIdx = 0;
    int searchCount = 0;
    while ( true ) {
        LinearNode & curNode = nodes[curNodeIdx];
        searchCount ++;
        if ( curNode.bounds.IntersectP(ray, invDir, dirIsNeg)   ) {
            // node case
            if ( curNode.nPrimitives > 0 ) {
                for ( int i = 0 ; i < curNode.nPrimitives ; ++ i ) {
                    auto its = primitives[curNode.primitivesOffset + i]->intersect(_ray);
                    if ( its.has_value() )
                        return true;
                }
                if ( visitIdx == 0 ) break; //no node to visit
                curNodeIdx = toVisit[-- visitIdx];
            }
            else {
                if ( dirIsNeg[curNode.axis] ) {
                    curNodeIdx += 1;
                    toVisit[visitIdx ++] = curNode.secondChildOffset;
                } else {
                    toVisit[visitIdx ++] = curNodeIdx + 1;
                    curNodeIdx = curNode.secondChildOffset;
                }
            }
        } else {
            if ( visitIdx == 0 ) break;
            curNodeIdx = toVisit[-- visitIdx];
        }
    }
    return false;
}

std::unique_ptr < BVHAccel > CreateBVH(const Json & j,
                                       const std::vector < std::shared_ptr < Primitive>> & primitives
) {

    return std::make_unique < BVHAccel >(primitives);
}