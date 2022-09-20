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

#include "subdivpatch1cached.h"
#include "grid_soa_intersector1.h"
#include "grid_soa_intersector.h"
#include "grid_aos_intersector1.h"
#include "../common/ray.h"

namespace embree
{
  namespace isa
  {
    class SubdivPatch1CachedIntersector1
    {
    public:
      typedef SubdivPatch1Cached Primitive;
      typedef GridSOAIntersector1::Precalculations Precalculations;
      
      static __forceinline bool processLazyNode(Precalculations& pre, const Primitive* prim_i, Scene* scene, size_t& lazy_node)
      {
        Primitive* prim = (Primitive*) prim_i;
        if (pre.grid) SharedLazyTessellationCache::sharedLazyTessellationCache.unlock();
        GridSOA* grid = (GridSOA*) SharedLazyTessellationCache::lookup(prim->entry(),scene->commitCounterSubdiv,[&] () {
            auto alloc = [] (const size_t bytes) { return SharedLazyTessellationCache::sharedLazyTessellationCache.malloc(bytes); };
            return GridSOA::create(prim,scene,alloc);
          });
        //GridSOA* grid = (GridSOA*) prim->root_ref.data;
        //GridSOA* grid = (GridSOA*) prim;
        lazy_node = grid->root;
        pre.grid = grid;
        return false;
      }

      /*! Intersect a ray with the primitive. */
      static __forceinline void intersect(Precalculations& pre, Ray& ray, const RTCIntersectContext* context, const Primitive* prim, size_t ty, Scene* scene, const unsigned* geomID_to_instID, size_t& lazy_node) 
      {
        if (likely(ty == 0)) GridSOAIntersector1::intersect(pre,ray,context,prim,ty,scene,lazy_node);
        else                 processLazyNode(pre,prim,scene,lazy_node);
      }
      static __forceinline void intersect(Precalculations& pre, Ray& ray, const RTCIntersectContext* context, size_t ty0, const Primitive* prim, size_t ty, Scene* scene, const unsigned* geomID_to_instID, size_t& lazy_node) {
        intersect(pre,ray,context,prim,ty,scene,geomID_to_instID,lazy_node);
      }
      
      /*! Test if the ray is occluded by the primitive */
      static __forceinline bool occluded(Precalculations& pre, Ray& ray, const RTCIntersectContext* context, const Primitive* prim, size_t ty, Scene* scene, const unsigned* geomID_to_instID, size_t& lazy_node) 
      {
        if (likely(ty == 0)) return GridSOAIntersector1::occluded(pre,ray,context,prim,ty,scene,lazy_node);
        else                 return processLazyNode(pre,prim,scene,lazy_node);
      }
      static __forceinline bool occluded(Precalculations& pre, Ray& ray, const RTCIntersectContext* context, size_t ty0, const Primitive* prim, size_t ty, Scene* scene, const unsigned* geomID_to_instID, size_t& lazy_node) {
        return occluded(pre,ray,context,prim,ty,scene,geomID_to_instID,lazy_node);
      }
    };

    template <int K>
    struct SubdivPatch1CachedIntersectorK
    {
      typedef SubdivPatch1Cached Primitive;
      typedef typename GridSOAIntersectorK<K>::Precalculations Precalculations;
      
      static __forceinline bool processLazyNode(Precalculations& pre, const Primitive* prim_i, Scene* scene, size_t& lazy_node)
      {
        Primitive* prim = (Primitive*) prim_i;
        if (pre.grid) SharedLazyTessellationCache::sharedLazyTessellationCache.unlock();
        GridSOA* grid = (GridSOA*) SharedLazyTessellationCache::lookup(prim->entry(),scene->commitCounterSubdiv,[&] () {
            auto alloc = [] (const size_t bytes) { return SharedLazyTessellationCache::sharedLazyTessellationCache.malloc(bytes); };
            return GridSOA::create(prim,scene,alloc);
          });
        lazy_node = grid->root;
        pre.grid = grid;
        return false;
      }
      
      static __forceinline void intersect(const vbool<K>& valid, Precalculations& pre, RayK<K>& ray, const RTCIntersectContext* context, const Primitive* prim, size_t ty, Scene* scene, size_t& lazy_node)
      {
        if (likely(ty == 0)) GridSOAIntersectorK<K>::intersect(valid,pre,ray,context,prim,ty,scene,lazy_node);
        else                 processLazyNode(pre,prim,scene,lazy_node);
      }
      
      static __forceinline vbool<K> occluded(const vbool<K>& valid, Precalculations& pre, RayK<K>& ray, const RTCIntersectContext* context, const Primitive* prim, size_t ty, Scene* scene, size_t& lazy_node)
      {
        if (likely(ty == 0)) return GridSOAIntersectorK<K>::occluded(valid,pre,ray,context,prim,ty,scene,lazy_node);
        else                 return processLazyNode(pre,prim,scene,lazy_node);
      }
      
      static __forceinline void intersect(Precalculations& pre, RayK<K>& ray, size_t k, const RTCIntersectContext* context, const Primitive* prim, size_t ty, Scene* scene, size_t& lazy_node)
      {
        if (likely(ty == 0)) GridSOAIntersectorK<K>::intersect(pre,ray,k,context,prim,ty,scene,lazy_node);
        else                 processLazyNode(pre,prim,scene,lazy_node);
      }
      
      static __forceinline bool occluded(Precalculations& pre, RayK<K>& ray, size_t k, const RTCIntersectContext* context, const Primitive* prim, size_t ty, Scene* scene, size_t& lazy_node)
      {
        if (likely(ty == 0)) return GridSOAIntersectorK<K>::occluded(pre,ray,k,context,prim,ty,scene,lazy_node);
        else                 return processLazyNode(pre,prim,scene,lazy_node);
      }
    };

    typedef SubdivPatch1CachedIntersectorK<4>  SubdivPatch1CachedIntersector4;
    typedef SubdivPatch1CachedIntersectorK<8>  SubdivPatch1CachedIntersector8;
    typedef SubdivPatch1CachedIntersectorK<16> SubdivPatch1CachedIntersector16;
  }
}
