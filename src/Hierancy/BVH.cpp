#include "BVH.hpp"


static int BuildNodeCount=0;

struct BuildNode
{
    BuildNode() { }


    void initLeaf(int off,int n,Bounds3 bounds3){
        firstPrimOffset = off;
        nPrimitives=n;
        bounds3=std::move(bounds3);
    }
    void    initInterior(int dim,BuildNode *  left, BuildNode * right )
    {
        nPrimitives=0;
        splitAxis=dim;
        children[0]=left;
        children[1]=right;
    }

    Bounds3 BB;
    BuildNode *children[2];
    std::vector<std::shared_ptr<Primitive>> primitives;
    int firstPrimOffset;
    int nPrimitives;
    int splitAxis;
};

struct BVHPrimitiveInfo{
    BVHPrimitiveInfo() {}
    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3 &bounds)
            : primitiveNumber(primitiveNumber),
              bounds(bounds),
              centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}
    size_t primitiveNumber;
    Bounds3 bounds;
    vec3 centroid;
};

struct alignas(32) LinearNode{
    Bounds3 bounds;
    union {
        int primitivesOffset;   // leaf
        int secondChildOffset;  // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;
};


BVHAccel::BVHAccel(const std::vector < std::shared_ptr < Primitive>> & primitives, const SplitMethod splitMethod,
                   const uint32 maxPrimInNode) {
    if(primitives.empty()) return ;

    //Take the primitive and its spatial information apart
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        primitiveInfo[i] = {i, primitives[i]->BB()};

    std::vector<std::shared_ptr<Primitive>> orderedPrims;

    BuildNode * node = RecursiveBuild(primitiveInfo,0,primitives.size(),orderedPrims);

    nodes.resize(BuildNodeCount);
    int offset=0;
    FlattenTree(node,&offset);
}

BuildNode  *  BVHAccel::RecursiveBuild( std::vector<BVHPrimitiveInfo> &primitiveInfo, int start, int end,
                            std::vector<std::shared_ptr<Primitive>> &orderedPrims ){
    if(start == end){
        return nullptr;
    }
    BuildNodeCount++;
    Bounds3 bounds;
    for(int i=start;i<end;i++){
        bounds = Union(bounds,primitiveInfo[i].bounds);
    }

    BuildNode * node = new BuildNode;
    int nPrimitives = end-start;

    //Only single primitive
    if(nPrimitives==1){
        int primIdx= orderedPrims.size();
        orderedPrims.push_back(primitives[primIdx]);
        node->initLeaf(primIdx,nPrimitives,bounds);
        return node;
    }

    Bounds3 centroidBounds;
    for (int i = start; i < end; ++i)
        centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
    int dim = centroidBounds.MaximumExtent();


    //All primitives are in one dimension
    if(centroidBounds.pMax[dim] == centroidBounds.pMin[dim]){
        int primIdx=primitiveInfo[start].primitiveNumber;
        orderedPrims.insert(orderedPrims.end(),&primitives[primIdx],&primitives[primIdx]+nPrimitives);
        node->initLeaf(primIdx,nPrimitives,bounds);
        return node;
    }

    int mid = (start+end)/2;

    switch ( splitMethod ) {
        case Middle:{
            Float pmid =(centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
            BVHPrimitiveInfo * midPtr = std::partition(
                    &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                    [dim, pmid](const BVHPrimitiveInfo &pi) {
                        return pi.centroid[dim] < pmid;
                    });
            mid = midPtr - &primitiveInfo[0];
            break;
        }
        case EqualCounts:{
            mid = (start + end) / 2;
            std::nth_element( &primitiveInfo[start], &primitiveInfo[mid],
                              &primitiveInfo[end - 1] + 1,
                              [dim](const BVHPrimitiveInfo &a,
                                    const BVHPrimitiveInfo &b) {
                                  return a.centroid[dim] < b.centroid[dim];
                              });}
            break;
        default:
            ; //todo
    }

    int primIdx= orderedPrims.size();
    node->initInterior(dim,
                       RecursiveBuild(primitiveInfo,start,mid,orderedPrims),
                       RecursiveBuild(primitiveInfo,mid,end,orderedPrims)
                       );

    return node;
}

//To be honest, I don't know why PBRT flatten the tree, it's probably for performance reasons,
// and I'll follow along
int  BVHAccel::FlattenTree(BuildNode * node, int * offset) {
    LinearNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->BB;
    (*offset)++;
    int myOffset=*offset;
    if(node->nPrimitives>0){
        linearNode->primitivesOffset=node->firstPrimOffset;
        linearNode->nPrimitives=node->nPrimitives;
    }
    else {
        linearNode->axis=node->splitAxis;
        FlattenTree(node->children[0],offset);
        linearNode->secondChildOffset = FlattenTree(node->children[1],offset);
    }
    return myOffset;

}

std::optional < Intersection > BVHAccel::intersect(const Ray & ray) {
    if(nodes.empty()){
        return std::nullopt;
    }

    Ray _ray(ray);

    vec3 invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};

    int toVisit[64];
    int visitIdx=0;
    int curNodeIdx=0;

    Intersection res;
    while(true){
        LinearNode & curNode=nodes[curNodeIdx];
        if(curNode.bounds.IntersectP(ray,invDir,dirIsNeg)){
            //left node case
            if(curNode.nPrimitives>0){
                for (int i = 0; i < curNode.nPrimitives; ++i)
                {
                    auto its=primitives[curNode.primitivesOffset + i]->intersect(_ray);
                    if(its.has_value()) res=its.value();
                }
                if(visitIdx==0) break; //no node to visit
            }
            //interior node case
            else {
                //near
                if(dirIsNeg[curNode.axis]){
                    curNodeIdx +=1;
                    toVisit[visitIdx++] = curNode.secondChildOffset;
                }
                else{
                    toVisit[visitIdx++] =curNodeIdx+1;
                    curNodeIdx=curNode.secondChildOffset;
                }
            }
        }
        else {
            curNodeIdx=toVisit[visitIdx--];
        }

    }




}

bool BVHAccel::intersectP(const Ray & ray) {

}

std::shared_ptr<BVHAccel> CreateBVH(const nlohmann::json &j,
                                    const std::vector<std::shared_ptr<Primitive>> primitives
){

    return std::make_shared <BVHAccel>(primitives);
}