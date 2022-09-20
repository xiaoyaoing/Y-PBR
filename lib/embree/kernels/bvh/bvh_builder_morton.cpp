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

#include "bvh.h"
#include "bvh_statistics.h"
#include "bvh_rotate.h"
#include "../common/profile.h"
#include "../algorithms/parallel_prefix_sum.h"
#include "../algorithms/sort.h"

#include "../builders/primrefgen.h"
#include "../builders/bvh_builder_morton.h"

#include "../geometry/triangle.h"
#include "../geometry/trianglev.h"
#include "../geometry/trianglei.h"

#define ROTATE_TREE 1 // specifies number of tree rotation rounds to perform

#define BLOCK_SIZE 1024

namespace embree 
{
  namespace isa
  {
    template<int N>
    struct AllocBVHNNode
    {
      typedef BVHN<N> BVH;
      typedef typename BVH::Node Node;
      typedef typename BVH::NodeRef NodeRef;

      __forceinline Node* operator() (MortonBuildRecord<NodeRef>& current, MortonBuildRecord<NodeRef>* children, size_t numChildren, FastAllocator::ThreadLocal2* alloc)
      {
        Node* node = (Node*) alloc->alloc0.malloc(sizeof(Node),BVH::byteNodeAlignment); 
        *current.parent = BVH::encodeNode(node);
        node->clear();
        for (size_t i=0; i<numChildren; i++)
          children[i].parent = &node->child(i);
        return node;
      }
    };
    
    template<int N>
    struct SetBVHNBounds
    {
      typedef BVHN<N> BVH;
      typedef typename BVH::Node Node;

      BVH* bvh;
      __forceinline SetBVHNBounds (BVH* bvh) : bvh(bvh) {}

      __forceinline BBox3fa operator() (Node* node, const BBox3fa* bounds, size_t num)
      {
        BBox3fa res = empty;
        for (size_t i=0; i<num; i++) {
          const BBox3fa b = bounds[i];
          res.extend(b);
          node->set(i,b);
        }

#if ROTATE_TREE
        if (N == 4)
        {
          size_t n = 0;
          for (size_t i=0; i<num; i++)
            n += bounds[i].lower.a;

          if (n >= 4096) {
            for (size_t i=0; i<num; i++) {
              if (bounds[i].lower.a < 4096) {
                for (int j=0; j<ROTATE_TREE; j++)
                  BVHNRotate<N>::rotate(node->child(i));
                node->child(i).setBarrier();
              }
            }
          }
          res.lower.a = unsigned(n);
        }
#endif

        return res;
      }
    };

    template<int N, typename Primitive>
    struct CreateMortonLeaf;

    template<int N>
    struct CreateMortonLeaf<N,Triangle4>
    {
      typedef BVHN<N> BVH;
      typedef typename BVH::NodeRef NodeRef;

      __forceinline CreateMortonLeaf (TriangleMesh* mesh, MortonID32Bit* morton)
        : mesh(mesh), morton(morton) {}

      __noinline void operator() (MortonBuildRecord<NodeRef>& current, FastAllocator::ThreadLocal2* alloc, BBox3fa& box_o)
      {
        vfloat4 lower(pos_inf);
        vfloat4 upper(neg_inf);
        size_t items = current.size();
        size_t start = current.begin;
        assert(items<=4);
        
        /* allocate leaf node */
        Triangle4* accel = (Triangle4*) alloc->alloc1.malloc(sizeof(Triangle4));
        *current.parent = BVH::encodeLeaf((char*)accel,1);
        vint4 vgeomID = -1, vprimID = -1;
        Vec3vf4 v0 = zero, v1 = zero, v2 = zero;

        for (size_t i=0; i<items; i++)
        {
          const unsigned index = morton[start+i].index;
          const unsigned primID = index; 
          const unsigned geomID = this->mesh->id;
          const TriangleMesh* mesh = this->mesh;
          const TriangleMesh::Triangle& tri = mesh->triangle(primID);
          const Vec3fa& p0 = mesh->vertex(tri.v[0]);
          const Vec3fa& p1 = mesh->vertex(tri.v[1]);
          const Vec3fa& p2 = mesh->vertex(tri.v[2]);
          lower = min(lower,(vfloat4)p0,(vfloat4)p1,(vfloat4)p2);
          upper = max(upper,(vfloat4)p0,(vfloat4)p1,(vfloat4)p2);
          vgeomID [i] = geomID;
          vprimID [i] = primID;
          v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
          v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
          v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
        }

        Triangle4::store_nt(accel,Triangle4(v0,v1,v2,vgeomID,vprimID));
        box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
#if ROTATE_TREE
        if (N == 4)
          box_o.lower.a = unsigned(current.size());
#endif
      }
    
    private:
      TriangleMesh* mesh;
      MortonID32Bit* morton;
    };
    
    template<int N>
    struct CreateMortonLeaf<N,Triangle4v>
    {
      typedef BVHN<N> BVH;
      typedef typename BVH::NodeRef NodeRef;

      __forceinline CreateMortonLeaf (TriangleMesh* mesh, MortonID32Bit* morton)
        : mesh(mesh), morton(morton) {}
      
      __noinline void operator() (MortonBuildRecord<NodeRef>& current, FastAllocator::ThreadLocal2* alloc, BBox3fa& box_o)
      {
        vfloat4 lower(pos_inf);
        vfloat4 upper(neg_inf);
        size_t items = current.size();
        size_t start = current.begin;
        assert(items<=4);
        
        /* allocate leaf node */
        Triangle4v* accel = (Triangle4v*) alloc->alloc1.malloc(sizeof(Triangle4v));
        *current.parent = BVH::encodeLeaf((char*)accel,1);
        
        vint4 vgeomID = -1, vprimID = -1;
        Vec3vf4 v0 = zero, v1 = zero, v2 = zero;

        for (size_t i=0; i<items; i++)
        {
          const unsigned index = morton[start+i].index;
          const unsigned primID = index; 
          const unsigned geomID = this->mesh->id;
          const TriangleMesh* mesh = this->mesh;
          const TriangleMesh::Triangle& tri = mesh->triangle(primID);
          const Vec3fa& p0 = mesh->vertex(tri.v[0]);
          const Vec3fa& p1 = mesh->vertex(tri.v[1]);
          const Vec3fa& p2 = mesh->vertex(tri.v[2]);
          lower = min(lower,(vfloat4)p0,(vfloat4)p1,(vfloat4)p2);
          upper = max(upper,(vfloat4)p0,(vfloat4)p1,(vfloat4)p2);
          vgeomID [i] = geomID;
          vprimID [i] = primID;
          v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
          v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
          v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
        }
        Triangle4v::store_nt(accel,Triangle4v(v0,v1,v2,vgeomID,vprimID));
        box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
#if ROTATE_TREE
        if (N == 4)
          box_o.lower.a = current.size();
#endif
      }
    private:
      TriangleMesh* mesh;
      MortonID32Bit* morton;
    };

    template<int N>
    struct CreateMortonLeaf<N,Triangle4i>
    {
      typedef BVHN<N> BVH;
      typedef typename BVH::NodeRef NodeRef;

      __forceinline CreateMortonLeaf (TriangleMesh* mesh, MortonID32Bit* morton)
        : mesh(mesh), morton(morton) {}
      
      __noinline void operator() (MortonBuildRecord<NodeRef>& current, FastAllocator::ThreadLocal2* alloc, BBox3fa& box_o)
      {
        vfloat4 lower(pos_inf);
        vfloat4 upper(neg_inf);
        size_t items = current.size();
        size_t start = current.begin;
        assert(items<=4);
        
        /* allocate leaf node */
        Triangle4i* accel = (Triangle4i*) alloc->alloc1.malloc(sizeof(Triangle4i));
        *current.parent = BVH::encodeLeaf((char*)accel,1);
        
        vint4 vgeomID = -1, vprimID = -1;
        Vec3f* v0[4] = { nullptr, nullptr, nullptr, nullptr };
        vint4 v1 = zero, v2 = zero;
        
        for (size_t i=0; i<items; i++)
        {
          const unsigned index = morton[start+i].index;
          const unsigned primID = index; 
          const unsigned geomID = this->mesh->id;
          const TriangleMesh* mesh = this->mesh;
          const TriangleMesh::Triangle& tri = mesh->triangle(primID);
          const Vec3fa& p0 = mesh->vertex(tri.v[0]);
          const Vec3fa& p1 = mesh->vertex(tri.v[1]);
          const Vec3fa& p2 = mesh->vertex(tri.v[2]);
          lower = min(lower,(vfloat4)p0,(vfloat4)p1,(vfloat4)p2);
          upper = max(upper,(vfloat4)p0,(vfloat4)p1,(vfloat4)p2);
          vgeomID[i] = geomID;
          vprimID[i] = primID;
          v0[i] = (Vec3f*) mesh->vertexPtr(tri.v[0]); 
          v1[i] = int(ssize_t((int*)   mesh->vertexPtr(tri.v[1])-(int*)v0[i])); 
          v2[i] = int(ssize_t((int*)   mesh->vertexPtr(tri.v[2])-(int*)v0[i])); 
        }
        
        for (size_t i=items; i<4; i++)
        {
          vgeomID[i] = -1;
          vprimID[i] = -1;
          v0[i] = v0[0];
          v1[i] = 0; 
          v2[i] = 0;
        }
        
        new (accel) Triangle4i(v0,v1,v2,vgeomID,vprimID);
        box_o = BBox3fa((Vec3fa)lower,(Vec3fa)upper);
#if ROTATE_TREE
        if (N == 4)
          box_o.lower.a = current.size();
#endif
      }
    private:
      TriangleMesh* mesh;
      MortonID32Bit* morton;
    };
    
    template<typename Mesh>
    struct CalculateMeshBounds
    {
      __forceinline CalculateMeshBounds (Mesh* mesh)
        : mesh(mesh) {}
      
      __forceinline const BBox3fa operator() (const MortonID32Bit& morton) {
        return mesh->bounds(morton.index);
      }
      
    private:
      Mesh* mesh;
    };        
    
    template<int N, typename Mesh, typename Primitive>
    class BVHNMeshBuilderMorton : public Builder
    {
      typedef BVHN<N> BVH;
      typedef typename BVH::Node Node;
      typedef typename BVH::NodeRef NodeRef;

    public:
      
      BVHNMeshBuilderMorton (BVH* bvh, Mesh* mesh, const size_t minLeafSize, const size_t maxLeafSize)
        : bvh(bvh), mesh(mesh), minLeafSize(minLeafSize), maxLeafSize(maxLeafSize), numPrimitives(0), morton(bvh->device) {}
      
      /*! Destruction */
      ~BVHNMeshBuilderMorton () {
        //bvh->shrink();
      }
      
      /* build function */
      void build(size_t threadIndex, size_t threadCount) 
      {
        /* We have to clear the allocator to guarantee that we can
         * temporarily use the first allocation block for sorting the
         * morton codes. */

        const size_t numNewPrimitives = mesh->size();
        if (numNewPrimitives != numPrimitives) bvh->alloc.clear();
        numPrimitives = numNewPrimitives;
        
        /* skip build for empty scene */
        if (numPrimitives == 0) {
          bvh->set(BVH::emptyNode,empty,0);
          return;
        }
        
        auto progress = [&] (size_t dn) { bvh->scene->progressMonitor(double(dn)); };
        
        /* preallocate arrays */
        morton.resize(numPrimitives);
        size_t bytesAllocated = numPrimitives*sizeof(Node)/(4*N) + size_t(1.2f*Primitive::blocks(numPrimitives)*sizeof(Primitive));
        size_t bytesMortonCodes = numPrimitives*sizeof(MortonID32Bit);
        bytesAllocated = max(bytesAllocated,bytesMortonCodes); // the first allocation block is reused to sort the morton codes
        bvh->alloc.init(bytesAllocated,2*bytesAllocated);

        size_t block_size = size_t(BLOCK_SIZE);

        /* compute scene bounds */
        std::pair<size_t,BBox3fa> cb_empty(0,empty);
        auto cb = parallel_reduce 
          ( size_t(0), numPrimitives, block_size, cb_empty, [&](const range<size_t>& r) -> std::pair<size_t,BBox3fa>
            {
              size_t num = 0;
              BBox3fa bounds = empty;

              for (size_t j=r.begin(); j<r.end(); j++)
              {
                BBox3fa prim_bounds = empty;
                if (unlikely(!mesh->valid(j,&prim_bounds))) continue;
                bounds.extend(center2(prim_bounds));
                num++;
              }
              return std::make_pair(num,bounds);
            }, [] (const std::pair<size_t,BBox3fa>& a, const std::pair<size_t,BBox3fa>& b) {
              return std::make_pair(a.first + b.first,merge(a.second,b.second)); 
            });

        size_t numPrimitivesGen = cb.first;
        const BBox3fa centBounds = cb.second;

        /* compute morton codes */
        MortonID32Bit* dest = (MortonID32Bit*) bvh->alloc.specialAlloc(bytesMortonCodes);

        if (likely(numPrimitivesGen == numPrimitives))
        {
          /* fast path */
          MortonCodeGenerator::MortonCodeMapping mapping(centBounds);
          parallel_for( size_t(0), numPrimitives, block_size, [&](const range<size_t>& r) -> void {
              MortonCodeGenerator generator(mapping,&morton.data()[r.begin()]);
              for (size_t j=r.begin(); j<r.end(); j++)
                generator(mesh->bounds(j),unsigned(j));
            });
        }
        else
        {
          /* slow path, fallback in case some primitives were invalid */
          ParallelPrefixSumState<size_t> pstate;
          MortonCodeGenerator::MortonCodeMapping mapping(centBounds);
          parallel_prefix_sum( pstate, size_t(0), numPrimitives, block_size, size_t(0), [&](const range<size_t>& r, const size_t base) -> size_t {
              size_t num = 0;
              MortonCodeGenerator generator(mapping,&morton.data()[r.begin()]);
              for (size_t j=r.begin(); j<r.end(); j++)
              {
                BBox3fa bounds = empty;
                if (unlikely(!mesh->valid(j,&bounds))) continue;
                generator(bounds,unsigned(j));
                num++;
              }
              return num;
            }, std::plus<size_t>());

          parallel_prefix_sum( pstate, size_t(0), numPrimitives, block_size, size_t(0), [&](const range<size_t>& r, const size_t base) -> size_t {
              size_t num = 0;
              MortonCodeGenerator generator(mapping,&morton.data()[base]);
              for (size_t j=r.begin(); j<r.end(); j++)
              {
                BBox3fa bounds = empty;
                if (!mesh->valid(j,&bounds)) continue;
                generator(bounds,unsigned(j));
                num++;
              }
              return num;
            }, std::plus<size_t>());          
        }

        /* create BVH */
        AllocBVHNNode<N> allocNode;
        SetBVHNBounds<N> setBounds(bvh);
        CreateMortonLeaf<N,Primitive> createLeaf(mesh,morton.data());
        CalculateMeshBounds<Mesh> calculateBounds(mesh);
        auto node_bounds = bvh_builder_morton_internal<NodeRef>(
          typename BVH::CreateAlloc(bvh), BBox3fa(empty),
          allocNode,setBounds,createLeaf,calculateBounds,progress,
          morton.data(),dest,numPrimitivesGen,N,BVH::maxBuildDepth,minLeafSize,maxLeafSize);
        
        bvh->set(node_bounds.first,node_bounds.second,numPrimitives);
        
#if ROTATE_TREE
        if (N == 4)
        {
          for (int i=0; i<ROTATE_TREE; i++)
            BVHNRotate<N>::rotate(bvh->root);
          bvh->clearBarrier(bvh->root);
        }
#endif

        /* clear temporary data for static geometry */
        if (mesh->isStatic()) 
        {
          morton.clear();
          bvh->shrink();
        }
        bvh->cleanup();
      }
      
      void clear() {
        morton.clear();
      }
      
    private:
      BVH* bvh;
      Mesh* mesh;
      const size_t minLeafSize;
      const size_t maxLeafSize;
      size_t numPrimitives;
      mvector<MortonID32Bit> morton;
    };

#if defined(EMBREE_GEOMETRY_TRIANGLES)
    Builder* BVH4Triangle4MeshBuilderMortonGeneral  (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVHNMeshBuilderMorton<4,TriangleMesh,Triangle4> ((BVH4*)bvh,mesh,4,4*BVH4::maxLeafBlocks); }
    Builder* BVH4Triangle4vMeshBuilderMortonGeneral (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVHNMeshBuilderMorton<4,TriangleMesh,Triangle4v>((BVH4*)bvh,mesh,4,4*BVH4::maxLeafBlocks); }
    Builder* BVH4Triangle4iMeshBuilderMortonGeneral (void* bvh, TriangleMesh* mesh, size_t mode) { return new class BVHNMeshBuilderMorton<4,TriangleMesh,Triangle4i>((BVH4*)bvh,mesh,4,4*BVH4::maxLeafBlocks); }
#endif
  }
}

