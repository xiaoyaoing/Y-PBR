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

#pragma once

#include "bvh.h"

namespace embree
{
  template<int N>
  class BVHNStatistics
  {
    typedef BVHN<N> BVH;
    typedef typename BVH::Node AlignedNode;
    typedef typename BVH::UnalignedNode UnalignedNode;
    typedef typename BVH::NodeMB AlignedNodeMB;
    typedef typename BVH::UnalignedNodeMB UnalignedNodeMB;
    typedef typename BVH::TransformNode TransformNode;
    typedef typename BVH::QuantizedNode QuantizedNode;

    typedef typename BVH::NodeRef NodeRef;

  public:

    /* Constructor gathers statistics. */
    BVHNStatistics (BVH* bvh);

    /*! Convert statistics into a string */
    std::string str();

    float sah() const { return bvhSAH; }

    size_t bytesUsed() const;

  private:
    void statistics(NodeRef node, const float A, size_t& depth);

  private:
    BVH* bvh;
    float bvhSAH;                      //!< SAH cost.
    float leafSAH;
    size_t numAlignedNodes;            //!< Number of aligned internal nodes.
    size_t numUnalignedNodes;          //!< Number of unaligned internal nodes.
    size_t numAlignedNodesMB;          //!< Number of aligned internal nodes.
    size_t numUnalignedNodesMB;        //!< Number of unaligned internal nodes.
    size_t numTransformNodes;          //!< Number of transformation nodes;
    size_t numQuantizedNodes;          //!< Number of transformation nodes;
    size_t childrenAlignedNodes;       //!< Number of children of aligned nodes
    size_t childrenUnalignedNodes;     //!< Number of children of unaligned internal nodes.
    size_t childrenAlignedNodesMB;     //!< Number of children of aligned nodes
    size_t childrenUnalignedNodesMB;   //!< Number of children of unaligned internal nodes.
    size_t childrenQuantizedNodes;     //!< Number of children of quantized internal nodes.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrims;                   //!< Number of primitives.
    size_t numPrimBlocks;              //!< Number of primitive blocks.
    size_t depth;                      //!< Depth of the tree.
  };

  typedef BVHNStatistics<4> BVH4Statistics;
  typedef BVHNStatistics<8> BVH8Statistics;
}
