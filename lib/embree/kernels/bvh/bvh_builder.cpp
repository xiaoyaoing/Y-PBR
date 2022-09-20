// ======================================================================== //
// Copyright 2009-2016 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "bvh_builder.h"
#include "bvh_rotate.h"

#define ROTATE_TREE 0

namespace embree
{
  namespace isa
  {
    /* tree rotations */
    template<int N>
    __forceinline size_t rotate(typename BVHN<N>::Node* node, const size_t* counts, const size_t num) {
      return 0;
    }

    template<int N>
    __forceinline size_t dummy(typename BVHN<N>::Node* node, const size_t* counts, const size_t num) {
      return 0;
    }

#if ROTATE_TREE
    template<>
    __forceinline size_t rotate<4>(BVH4::Node* node, const size_t* counts, const size_t num)
    {
      size_t n = 0;
      assert(num <= 4);
      for (size_t i=0; i<num; i++)
        n += counts[i];
      if (n >= 4096) {
        for (size_t i=0; i<num; i++) {
          if (counts[i] < 4096) {
            for (int j=0; j<ROTATE_TREE; j++) 
              BVHNRotate<4>::rotate(node->child(i));
            node->child(i).setBarrier();
          }
        }
      }
      return n;
    }
#endif

    template<int N>
    void BVHNBuilder<N>::BVHNBuilderV::build(BVH* bvh, BuildProgressMonitor& progress_in, PrimRef* prims, const PrimInfo& pinfo, const size_t blockSize, const size_t minLeafSize, const size_t maxLeafSize, const float travCost, const float intCost)
    {
      //bvh->alloc.init_estimate(pinfo.size()*sizeof(PrimRef));

      auto progressFunc = [&] (size_t dn) { 
        progress_in(dn); 
      };
            
      auto createLeafFunc = [&] (const BVHBuilderBinnedSAH::BuildRecord& current, Allocator* alloc) -> size_t {
        return createLeaf(current,alloc);
      };
      
      NodeRef root;
      BVHBuilderBinnedSAH::build_reduce<NodeRef>
        (root,typename BVH::CreateAlloc(bvh),size_t(0),typename BVH::CreateNode(bvh),rotate<N>,createLeafFunc,progressFunc,
         prims,pinfo,N,BVH::maxBuildDepthLeaf,blockSize,minLeafSize,maxLeafSize,travCost,intCost);

      bvh->set(root,pinfo.geomBounds,pinfo.size());
      
#if ROTATE_TREE
      if (N == 4)
      {
        for (int i=0; i<ROTATE_TREE; i++)
          BVHNRotate<N>::rotate(bvh->root);
        bvh->clearBarrier(bvh->root);
      }
#endif
      
      bvh->layoutLargeNodes(size_t(pinfo.size()*0.005f));
    }


    template<int N>
    void BVHNBuilderQuantized<N>::BVHNBuilderV::build(BVH* bvh, BuildProgressMonitor& progress_in, PrimRef* prims, const PrimInfo& pinfo, const size_t blockSize, const size_t minLeafSize, const size_t maxLeafSize, const float travCost, const float intCost)
    {
      //bvh->alloc.init_estimate(pinfo.size()*sizeof(PrimRef));
      auto progressFunc = [&] (size_t dn) { 
        progress_in(dn); 
      };
            
      auto createLeafFunc = [&] (const BVHBuilderBinnedSAH::BuildRecord& current, Allocator* alloc) -> size_t {
        return createLeaf(current,alloc);
      };
            
#if ENABLE_32BIT_OFFSETS_FOR_QUANTIZED_NODES == 1 
      typename BVH::QuantizedNode *first = (typename BVH::QuantizedNode*)bvh->alloc.malloc(sizeof(typename BVH::QuantizedNode), bvh->byteNodeAlignment);
      NodeRef &root = *(NodeRef*)first; // as the builder assigns current.parent = (size_t*)&root
      assert(((size_t)first & 0x7) == 0); 
#else
      NodeRef root = 0;
#endif
      BVHBuilderBinnedSAH::build_reduce<NodeRef>
        (root,typename BVH::CreateAlloc(bvh),size_t(0),typename BVH::CreateQuantizedNode(bvh),dummy<N>,createLeafFunc,progressFunc,
         prims,pinfo,N,BVH::maxBuildDepthLeaf,blockSize,minLeafSize,maxLeafSize,travCost,intCost);

#if ENABLE_32BIT_OFFSETS_FOR_QUANTIZED_NODES == 1 
      NodeRef new_root = ((size_t)first + first->childOffset(0)) | BVH::tyQuantizedNode;
#else
      NodeRef new_root = (size_t)root | BVH::tyQuantizedNode;
      // todo: COPY LAYOUT FOR LARGE NODES !!!
      //bvh->layoutLargeNodes(pinfo.size()*0.005f);
#endif
      assert(new_root.isQuantizedNode());
      bvh->set(new_root,pinfo.geomBounds,pinfo.size());
    }

    template<int N>
      struct CreateNodeMB
    {
      typedef BVHN<N> BVH;
      typedef typename BVH::NodeMB NodeMB;

      __forceinline CreateNodeMB (BVH* bvh) : bvh(bvh) {}
      
      __forceinline NodeMB* operator() (const isa::BVHBuilderBinnedSAH::BuildRecord& current, BVHBuilderBinnedSAH::BuildRecord* children, const size_t num, FastAllocator::ThreadLocal2* alloc)
      {
        NodeMB* node = (NodeMB*) alloc->alloc0.malloc(sizeof(NodeMB)); node->clear();
        for (size_t i=0; i<num; i++) {
          children[i].parent = (size_t*)&node->child(i);
        }
        *current.parent = bvh->encodeNode(node);
	return node;
      }

      BVH* bvh;
    };

    template<int N>
    void BVHNBuilderMblur<N>::BVHNBuilderV::build(BVH* bvh, BuildProgressMonitor& progress_in, PrimRef* prims, const PrimInfo& pinfo, const size_t blockSize, const size_t minLeafSize, const size_t maxLeafSize, const float travCost, const float intCost)
    {
      //bvh->alloc.init_estimate(pinfo.size()*sizeof(PrimRef));

      auto progressFunc = [&] (size_t dn) { 
        progress_in(dn); 
      };
            
      auto createLeafFunc = [&] (const BVHBuilderBinnedSAH::BuildRecord& current, Allocator* alloc) -> std::pair<BBox3fa,BBox3fa> {
        return createLeaf(current,alloc);
      };

      /* reduction function */
      auto reduce = [] (NodeMB* node, const std::pair<BBox3fa,BBox3fa>* bounds, const size_t num) -> std::pair<BBox3fa,BBox3fa>
      {
        assert(num <= N);
        BBox3fa bounds0 = empty;
        BBox3fa bounds1 = empty;
        for (size_t i=0; i<num; i++) {
          const BBox3fa b0 = bounds[i].first;
          const BBox3fa b1 = bounds[i].second;
          node->set(i,b0,b1);
          bounds0 = merge(bounds0,b0);
          bounds1 = merge(bounds1,b1);
        }
        return std::pair<BBox3fa,BBox3fa>(bounds0,bounds1);
      };
      auto identity = std::make_pair(BBox3fa(empty),BBox3fa(empty));
      
      
      NodeRef root;
      std::pair<BBox3fa,BBox3fa> root_bounds = BVHBuilderBinnedSAH::build_reduce<NodeRef>
        (root,typename BVH::CreateAlloc(bvh),identity,CreateNodeMB<N>(bvh),reduce,createLeafFunc,progressFunc,
         prims,pinfo,N,BVH::maxBuildDepthLeaf,blockSize,minLeafSize,maxLeafSize,travCost,intCost);

      /* set bounding box to merge bounds of all time steps */
      //bvh->set(root,pinfo.geomBounds,pinfo.size());
      bvh->set(root,merge(root_bounds.first,root_bounds.second),pinfo.size());


#if ROTATE_TREE
      if (N == 4)
      {
        for (int i=0; i<ROTATE_TREE; i++)
          BVHNRotate<N>::rotate(bvh->root);
        bvh->clearBarrier(bvh->root);
      }
#endif
      
      //bvh->layoutLargeNodes(pinfo.size()*0.005f); // FIXME: implement for Mblur nodes and activate
    }

    template<int N>
    void BVHNBuilderSpatial<N>::BVHNBuilderV::build(BVH* bvh, BuildProgressMonitor& progress_in, 
                                                    PrimRefList& prims, const PrimInfo& pinfo, 
                                                    const size_t blockSize, const size_t minLeafSize, 
                                                    const size_t maxLeafSize, const float travCost, 
                                                    const float intCost)
    {
      //bvh->alloc.init_estimate(pinfo.size()*sizeof(PrimRef));
      
      auto progressFunc = [&] (size_t dn) { 
        progress_in(dn); 
      };

      auto splitPrimitiveFunc = [&] (const PrimRef& prim, int dim, float pos, PrimRef& left_o, PrimRef& right_o) -> void {
        splitPrimitive(prim,dim,pos,left_o,right_o);
      };

      auto createLeafFunc = [&] (BVHBuilderBinnedSpatialSAH::BuildRecord& current, Allocator* alloc) -> size_t {
        return createLeaf(current,alloc);
      };
      
      NodeRef root;


      BVHBuilderBinnedSpatialSAH::build_reduce<NodeRef>
        (root,typename BVH::CreateAlloc(bvh),size_t(0),typename BVH::CreateNode(bvh),rotate<N>,
         createLeafFunc,splitPrimitiveFunc,progressFunc,
         prims,pinfo,N,BVH::maxBuildDepthLeaf,blockSize,minLeafSize,maxLeafSize,travCost,intCost);
      
      bvh->set(root,pinfo.geomBounds,pinfo.size());
      
#if ROTATE_TREE
      if (N == 4)
      {
        for (int i=0; i<ROTATE_TREE; i++)
          BVHNRotate<N>::rotate(bvh->root);
        bvh->clearBarrier(bvh->root);
      }
#endif
      
      bvh->layoutLargeNodes(size_t(pinfo.size()*0.005f));
    }

    // ========================================================================================================================================================
    // ========================================================================================================================================================
    // ========================================================================================================================================================

    template<int N, int NUM_SPATIAL_SPLITS>
    void BVHNBuilderFastSpatial<N,NUM_SPATIAL_SPLITS>::BVHNBuilderV::build(BVH* bvh, 
                                                                           BuildProgressMonitor& progress_in, 
                                                                           PrimRef* prims0, 
                                                                           const size_t extSize,
                                                                           const PrimInfo& pinfo, const size_t blockSize, 
                                                                           const size_t minLeafSize, const size_t maxLeafSize, 
                                                                           const float travCost, const float intCost)
    {
      auto progressFunc = [&] (size_t dn) { 
        progress_in(dn); 
      };

      auto splitPrimitiveFunc = [&] (const PrimRef& prim, int dim, float pos, PrimRef& left_o, PrimRef& right_o) -> void {
        splitPrimitive(prim,dim,pos,left_o,right_o);
      };

      auto binnerSplitPrimitiveFunc = [&] (SpatialBinInfo<NUM_SPATIAL_SPLITS,PrimRef> &spatialBinner, const PrimRef* const source, const size_t begin, const size_t end, const SpatialBinMapping<NUM_SPATIAL_SPLITS> &mapping) -> void {
        binnerSplit(spatialBinner,source,begin,end,mapping);
      };
            
      auto createLeafFunc = [&] (const BVHBuilderBinnedFastSpatialSAH::BuildRecord& current, Allocator* alloc) -> size_t {
        return createLeaf(current,alloc);
      };
      
      NodeRef root;
      BVHBuilderBinnedFastSpatialSAH::build_reduce<NodeRef>
        (root,typename BVH::CreateAlloc(bvh),size_t(0),typename BVH::CreateNode(bvh),rotate<N>,createLeafFunc,splitPrimitiveFunc,binnerSplitPrimitiveFunc, progressFunc,
         prims0,extSize,pinfo,N,BVH::maxBuildDepthLeaf,blockSize,minLeafSize,maxLeafSize,travCost,intCost);

      bvh->set(root,pinfo.geomBounds,pinfo.size());      
      bvh->layoutLargeNodes(size_t(pinfo.size()*0.005f));
    }


    template struct BVHNBuilder<4>;
    template struct BVHNBuilderQuantized<4>;
    template struct BVHNBuilderMblur<4>;    
    template struct BVHNBuilderSpatial<4>;
    template struct BVHNBuilderFastSpatial<4,FAST_SPATIAL_BUILDER_NUM_SPATIAL_SPLITS>;

#if defined(__AVX__)
    template struct BVHNBuilder<8>;
    template struct BVHNBuilderQuantized<8>;
    template struct BVHNBuilderMblur<8>;
    template struct BVHNBuilderSpatial<8>;
    template struct BVHNBuilderFastSpatial<8,FAST_SPATIAL_BUILDER_NUM_SPATIAL_SPLITS>;
#endif
  }
}
