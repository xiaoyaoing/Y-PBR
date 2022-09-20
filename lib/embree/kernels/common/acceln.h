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

#include "accel.h"

namespace embree
{
  /*! merges N acceleration structures together, by processing them in order */
  class AccelN : public Accel
  {
  public:
    AccelN ();
    ~AccelN();

  public:
    void add(Accel* accel);

  public:
    static void intersect (void* ptr, RTCRay& ray, const RTCIntersectContext* context);
    static void intersect4 (const void* valid, void* ptr, RTCRay4& ray, const RTCIntersectContext* context);
    static void intersect8 (const void* valid, void* ptr, RTCRay8& ray, const RTCIntersectContext* context);
    static void intersect16 (const void* valid, void* ptr, RTCRay16& ray, const RTCIntersectContext* context);
    static void intersectN (void* ptr, RTCRay** ray, const size_t N, const RTCIntersectContext* context);

  public:
    static void occluded (void* ptr, RTCRay& ray, const RTCIntersectContext* context);
    static void occluded4 (const void* valid, void* ptr, RTCRay4& ray, const RTCIntersectContext* context);
    static void occluded8 (const void* valid, void* ptr, RTCRay8& ray, const RTCIntersectContext* context);
    static void occluded16 (const void* valid, void* ptr, RTCRay16& ray, const RTCIntersectContext* context);
    static void occludedN (void* ptr, RTCRay** ray, const size_t N, const RTCIntersectContext* context);

  public:
    void print(size_t ident);
    void immutable();
    void build (size_t threadIndex, size_t threadCount);
    void select(bool filter4, bool filter8, bool filter16, bool filterN);
    void deleteGeometry(size_t geomID);
    void clear ();
      
  public:
    darray_t<Accel*,16> accels;
    darray_t<Accel*,16> validAccels;
  };
}
