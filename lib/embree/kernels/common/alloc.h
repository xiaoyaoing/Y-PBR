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

#include "default.h"

/* only activate in 64bit mode */
#if defined(__X86_64__) && !defined(__WIN32__) // FIXME: no idea why this doesn't work on windows
#define ENABLE_PARALLEL_BLOCK_ALLOCATION 1
#else
#define ENABLE_PARALLEL_BLOCK_ALLOCATION 0
#endif

namespace embree
{
  class FastAllocator 
  {
    /*! maximal supported alignment */
    static const size_t maxAlignment = 64;

    /*! maximal allocation size */

    /* default settings */
    static const size_t defaultBlockSize = 4096;
    static const size_t maxAllocationSize = 4*1024*1024-maxAlignment;
    static const size_t MAX_THREAD_USED_BLOCK_SLOTS = 8;
    
  public:

    /*! Per thread structure holding the current memory block. */
    struct __aligned(64) ThreadLocal 
    {
      ALIGNED_CLASS_(64);
    public:

      /*! Constructor for usage with ThreadLocalData */
      __forceinline ThreadLocal (void* alloc) 
	: alloc((FastAllocator*)alloc), ptr(nullptr), cur(0), end(0), allocBlockSize(defaultBlockSize), bytesUsed(0), bytesWasted(0) {}

      /*! Default constructor. */
      __forceinline ThreadLocal (FastAllocator* alloc, const size_t allocBlockSize = defaultBlockSize) 
	: alloc(alloc), ptr(nullptr), cur(0), end(0), allocBlockSize(allocBlockSize), bytesUsed(0), bytesWasted(0)  {}

      /*! resets the allocator */
      __forceinline void reset() 
      {
	ptr = nullptr;
	cur = end = 0;
	bytesWasted = bytesUsed = 0;
      }

      /* Allocate aligned memory from the threads memory block. */
      __forceinline void* operator() (size_t bytes, size_t align = 16) {
        return malloc(bytes,align);
      }

      /* Allocate aligned memory from the threads memory block. */
      __forceinline void* malloc(size_t bytes, size_t align = 16) 
      {
        assert(align <= maxAlignment);
	bytesUsed += bytes;
	
        /* try to allocate in local block */
	size_t ofs = (align - cur) & (align-1); 
        cur += bytes + ofs;
        if (likely(cur <= end)) { bytesWasted += ofs; return &ptr[cur - bytes]; }
	cur -= bytes + ofs;

        /* if allocation is too large allocate with parent allocator */
        if (4*bytes > allocBlockSize) {
          return alloc->malloc(bytes,maxAlignment);
	}

        /* get new full block if allocation failed */
        size_t blockSize = allocBlockSize;
	ptr = (char*) alloc->malloc(blockSize,maxAlignment);
	bytesWasted += end-cur;
	cur = 0; end = blockSize;
	
        /* retry allocation */
	ofs = (align - cur) & (align-1); 
        cur += bytes + ofs;
        if (likely(cur <= end)) { bytesWasted += ofs; return &ptr[cur - bytes]; }
	cur -= bytes + ofs;
	
        /* should never happen as large allocations get handled specially above */
        assert(false);
        return nullptr;
      }

      /* returns current address */
      __forceinline void* curPtr() {
        if (ptr == nullptr) ptr = (char*) alloc->malloc(allocBlockSize,maxAlignment);
        return &ptr[bytesUsed];
      }

      /*! returns amount of used bytes */
      size_t getUsedBytes() const { return bytesUsed; }
      
      /*! returns amount of wasted bytes */
      size_t getWastedBytes() const { return bytesWasted + (end-cur); }

    public:
      FastAllocator* alloc;  //!< parent allocator
      char*  ptr;            //!< pointer to memory block
      size_t cur;            //!< current location of the allocator
      size_t end;            //!< end of the memory block
      size_t allocBlockSize; //!< block size for allocations
    private:
      size_t bytesUsed;      //!< number of total bytes allocated
      size_t bytesWasted;    //!< number of bytes wasted
    };

    /*! Two thread local structures. */
    struct __aligned(64) ThreadLocal2
    {
      ALIGNED_STRUCT;

      /*! Constructor for usage with ThreadLocalData */
      __forceinline ThreadLocal2 (void* alloc) 
        : alloc0(alloc), alloc1(alloc) {}

      /*! Default constructor. */
      __forceinline ThreadLocal2 (FastAllocator* alloc, const size_t allocBlockSize = defaultBlockSize) 
        : alloc0(alloc,allocBlockSize), alloc1(alloc,allocBlockSize) {}

      /*! resets the allocator */
      __forceinline void reset() {
        alloc0.reset();
        alloc1.reset();
      }

      /*! returns amount of used bytes */
      size_t getUsedBytes() const { return alloc0.getUsedBytes() + alloc1.getUsedBytes(); }
      
      /*! returns amount of wasted bytes */
      size_t getWastedBytes() const { return alloc0.getWastedBytes() + alloc1.getWastedBytes(); }
    
    public:  
      ThreadLocal alloc0;
      ThreadLocal alloc1;
    };

    FastAllocator (MemoryMonitorInterface* device) 
      : device(device), slotMask(0), usedBlocks(nullptr), freeBlocks(nullptr), growSize(defaultBlockSize), bytesUsed(0), thread_local_allocators(this), thread_local_allocators2(this)
    {
      for (size_t i=0; i<MAX_THREAD_USED_BLOCK_SLOTS; i++)
      {
        threadUsedBlocks[i] = nullptr;
        threadBlocks[i] = nullptr;
        assert(!slotMutex[i].isLocked());
      }
    }

    ~FastAllocator () { 
      clear();
    }

    /*! returns a fast thread local allocator */
    __forceinline ThreadLocal* threadLocal() {
      return thread_local_allocators.get();
    }

    /*! returns a fast thread local allocator */
    __forceinline ThreadLocal2* threadLocal2() {
      return thread_local_allocators2.get();
    }

    /*! initializes the allocator */
    void init(size_t bytesAllocate, size_t bytesReserve = 0) 
    {     
      /* distribute the allocation to multiple thread block slots */
      slotMask = MAX_THREAD_USED_BLOCK_SLOTS-1;      
      if (usedBlocks.load() || freeBlocks.load()) { reset(); return; }
      if (bytesReserve == 0) bytesReserve = bytesAllocate;
      freeBlocks = Block::create(device,bytesAllocate,bytesReserve);
      growSize = max(size_t(defaultBlockSize),bytesReserve);
    }

    /*! initializes the allocator */
    void init_estimate(size_t bytesAllocate) 
    {
      if (usedBlocks.load() || freeBlocks.load()) { reset(); return; }
      growSize = max(size_t(defaultBlockSize),bytesAllocate);
      if (bytesAllocate > 4*maxAllocationSize) slotMask = 0x1;
      if (bytesAllocate > 16*maxAllocationSize) slotMask = MAX_THREAD_USED_BLOCK_SLOTS-1;
    }

    /*! frees state not required after build */
    __forceinline void cleanup() 
    {

#if ENABLE_PARALLEL_BLOCK_ALLOCATION == 1
      /* move thread local blocks to global block list */
      for (size_t i=0;i<MAX_THREAD_USED_BLOCK_SLOTS;i++)
      {
        while (threadBlocks[i].load() != nullptr) {
          Block* nextUsedBlock = threadBlocks[i].load()->next;
          threadBlocks[i].load()->next = usedBlocks.load();
          usedBlocks = threadBlocks[i].load();
          threadBlocks[i] = nextUsedBlock;
        }
        threadBlocks[i] = nullptr;
      }
#endif
      
      for (size_t t=0; t<thread_local_allocators.threads.size(); t++)
        bytesUsed += thread_local_allocators.threads[t]->getUsedBytes();

      for (size_t t=0; t<thread_local_allocators2.threads.size(); t++)
	bytesUsed += thread_local_allocators2.threads[t]->getUsedBytes();

      thread_local_allocators.clear();
      thread_local_allocators2.clear();
    }

    /*! shrinks all memory blocks to the actually used size */
    void shrink () {
      if (usedBlocks.load() != nullptr) usedBlocks.load()->shrink(device);
      if (freeBlocks.load() != nullptr) freeBlocks.load()->clear(device); freeBlocks = nullptr;
    }


    /*! resets the allocator, memory blocks get reused */
    void reset () 
    {
      bytesUsed = 0;

      /* first reset all used blocks */
      if (usedBlocks.load() != nullptr) usedBlocks.load()->reset();
      
      /* move all used blocks to begin of free block list */
      while (usedBlocks.load() != nullptr) {
        Block* nextUsedBlock = usedBlocks.load()->next;
        usedBlocks.load()->next = freeBlocks.load();
        freeBlocks = usedBlocks.load();
        usedBlocks = nextUsedBlock;
      }

      for (size_t i=0; i<MAX_THREAD_USED_BLOCK_SLOTS; i++) 
      {
        threadUsedBlocks[i] = nullptr;
        threadBlocks[i] = nullptr;
      }
      
      /* reset all thread local allocators */
      thread_local_allocators.reset();
      thread_local_allocators2.reset();
    }

    /*! frees all allocated memory */
    __forceinline void clear()
    {
      cleanup();
      bytesUsed = 0;
      if (usedBlocks.load() != nullptr) usedBlocks.load()->clear(device); usedBlocks = nullptr;
      if (freeBlocks.load() != nullptr) freeBlocks.load()->clear(device); freeBlocks = nullptr;
      for (size_t i=0; i<MAX_THREAD_USED_BLOCK_SLOTS; i++) {
        threadUsedBlocks[i] = nullptr;
        threadBlocks[i] = nullptr;
      }
    }



    /*! thread safe allocation of memory */
    __noinline void* malloc(size_t bytes, size_t align) 
    {
      assert(align <= maxAlignment);

      while (true) 
      {
        /* allocate using current block */
        size_t threadIndex = TaskScheduler::threadIndex();
        size_t slot = threadIndex & slotMask;
	Block* myUsedBlocks = threadUsedBlocks[slot];
        if (myUsedBlocks) {
          void* ptr = myUsedBlocks->malloc(device,bytes,align); 
          if (ptr) return ptr;
        }
        

        /* throw error if allocation is too large */
        if (bytes > maxAllocationSize)
          throw_RTCError(RTC_UNKNOWN_ERROR,"allocation is too large");

#if ENABLE_PARALLEL_BLOCK_ALLOCATION == 1
        /* parallel block creation in case of no freeBlocks, avoids single global mutex */
        if (likely(freeBlocks.load() == nullptr)) 
          {
            {
              Lock<SpinLock> lock(slotMutex[slot]);
              if (myUsedBlocks == threadUsedBlocks[slot])
              {
                threadBlocks[slot] = threadUsedBlocks[slot] = Block::create(device,maxAllocationSize-maxAlignment, maxAllocationSize-maxAlignment, threadBlocks[slot]);              
              }    
            }
            continue;
          }        
#endif

        /* if this fails allocate new block */
        {
          Lock<SpinLock> lock(mutex);
	  if (myUsedBlocks == threadUsedBlocks[slot])
	  {
            if (freeBlocks.load() != nullptr) {
	      Block* nextFreeBlock = freeBlocks.load()->next;
	      freeBlocks.load()->next = usedBlocks;
	      __memory_barrier();
	      usedBlocks = freeBlocks.load();
              threadUsedBlocks[slot] = freeBlocks.load();
	      freeBlocks = nextFreeBlock;
	    } else {
	      growSize = min(2*growSize,size_t(maxAllocationSize+maxAlignment));
	      usedBlocks = threadUsedBlocks[slot] = Block::create(device,growSize-maxAlignment, growSize-maxAlignment, usedBlocks);
	    }
	  }
        }
      }
    }

    /* special allocation only used from morton builder only a single time for each build */
    void* specialAlloc(size_t bytes) 
    {
      /* create a new block if the first free block is too small */
      if (freeBlocks.load() == nullptr || freeBlocks.load()->getBlockAllocatedBytes() < bytes)
        freeBlocks = Block::create(device,bytes,bytes,freeBlocks);

      return freeBlocks.load()->ptr();
    }

    size_t getAllocatedBytes() const 
    {
      size_t bytesAllocated = 0;
      if (freeBlocks.load()) bytesAllocated += freeBlocks.load()->getAllocatedBytes();
      if (usedBlocks.load()) bytesAllocated += usedBlocks.load()->getAllocatedBytes();
      return bytesAllocated;
    }

    size_t getReservedBytes() const 
    {
      size_t bytesReserved = 0;
      if (freeBlocks.load()) bytesReserved += freeBlocks.load()->getReservedBytes();
      if (usedBlocks.load()) bytesReserved += usedBlocks.load()->getReservedBytes();
      return bytesReserved;
    }

    size_t getUsedBytes() const 
    {
      size_t bytes = bytesUsed;

      for (size_t t=0; t<thread_local_allocators.threads.size(); t++)
	bytes += thread_local_allocators.threads[t]->getUsedBytes();

      for (size_t t=0; t<thread_local_allocators2.threads.size(); t++)
	bytes += thread_local_allocators2.threads[t]->getUsedBytes();

      return bytes;
    }

    size_t getFreeBytes() const 
    {
      size_t bytesFree = 0;
      if (freeBlocks.load()) bytesFree += freeBlocks.load()->getAllocatedBytes();
      if (usedBlocks.load()) bytesFree += usedBlocks.load()->getFreeBytes();
      return bytesFree;
    }

    size_t getWastedBytes() const 
    {
      size_t bytesWasted = 0;
      if (usedBlocks.load()) {
	Block* cur = usedBlocks;
	while ((cur = cur->next) != nullptr)
	  bytesWasted += cur->getFreeBytes();
      }
      
      for (size_t t=0; t<thread_local_allocators.threads.size(); t++)
	bytesWasted += thread_local_allocators.threads[t]->getWastedBytes();

      for (size_t t=0; t<thread_local_allocators2.threads.size(); t++)
	bytesWasted += thread_local_allocators2.threads[t]->getWastedBytes();
      
      return bytesWasted;
    }

    void print_statistics()
    {
      size_t bytesFree = getFreeBytes();
      size_t bytesAllocated = getAllocatedBytes();
      size_t bytesReserved = getReservedBytes();
      size_t bytesUsed = getUsedBytes();
      size_t bytesWasted = getWastedBytes();
      printf("  allocated = %3.2fMB, reserved = %3.2fMB, used = %3.2fMB (%3.2f%%), wasted = %3.2fMB (%3.2f%%), free = %3.2fMB (%3.2f%%)\n",
	     1E-6f*bytesAllocated, 1E-6f*bytesReserved,
	     1E-6f*bytesUsed, 100.0f*bytesUsed/bytesAllocated,
	     1E-6f*bytesWasted, 100.0f*bytesWasted/bytesAllocated,
	     1E-6f*bytesFree, 100.0f*bytesFree/bytesAllocated);
      
      //if (State::instance()->verbosity(3)) 
      {
        std::cout << "  used blocks = ";
        if (usedBlocks.load() != nullptr) usedBlocks.load()->print();
        std::cout << "[END]" << std::endl;
        
        std::cout << "  free blocks = ";
        if (freeBlocks.load() != nullptr) freeBlocks.load()->print();
        std::cout << "[END]" << std::endl;
      }
    }

  private:

    struct Block 
    {
      static Block* create(MemoryMonitorInterface* device, size_t bytesAllocate, size_t bytesReserve, Block* next = nullptr)
      {
        const size_t sizeof_Header = offsetof(Block,data[0]);
        bytesAllocate = ((sizeof_Header+bytesAllocate+defaultBlockSize-1) & ~(defaultBlockSize-1)); // always consume full pages
        bytesReserve  = ((sizeof_Header+bytesReserve +defaultBlockSize-1) & ~(defaultBlockSize-1)); // always consume full pages
        if (device) device->memoryMonitor(bytesAllocate,false);
        void* ptr = os_reserve(bytesReserve);
        os_commit(ptr,bytesAllocate);
        new (ptr) Block(bytesAllocate-sizeof_Header,bytesReserve-sizeof_Header,next);
        return (Block*) ptr;
      }

      Block (size_t bytesAllocate, size_t bytesReserve, Block* next) 
      : cur(0), allocEnd(bytesAllocate), reserveEnd(bytesReserve), next(next) 
      {
        //for (size_t i=0; i<allocEnd; i+=defaultBlockSize) data[i] = 0;
      }

      void clear (MemoryMonitorInterface* device) {
	if (next) next->clear(device); next = nullptr;
        const size_t sizeof_Header = offsetof(Block,data[0]);
        size_t sizeof_This = sizeof_Header+reserveEnd;
        const ssize_t sizeof_Alloced = sizeof_Header+getBlockAllocatedBytes();
        os_free(this,sizeof_This);
        if (device) device->memoryMonitor(-sizeof_Alloced,true);
      }
      
      void* malloc(MemoryMonitorInterface* device, size_t bytes, size_t align = 16) 
      {
        assert(align <= maxAlignment);
        bytes = (bytes+(align-1)) & ~(align-1);
	if (unlikely(cur+bytes > reserveEnd)) return nullptr;
	const size_t i = cur.fetch_add(bytes);
	if (unlikely(i+bytes > reserveEnd)) return nullptr;
	if (i+bytes > allocEnd) {
          if (device) device->memoryMonitor(i+bytes-max(i,allocEnd),true);
          os_commit(&data[i],bytes); // FIXME: optimize, may get called frequently
        }
	return &data[i];
      }
      
      void* ptr() {
        return &data[cur];
      }

      void reset () 
      {
        allocEnd = max(allocEnd,(size_t)cur);
        cur = 0;
        if (next) next->reset(); 
      }

      void shrink (MemoryMonitorInterface* device) 
      {
        const size_t sizeof_Header = offsetof(Block,data[0]);
        size_t newSize = os_shrink(this,sizeof_Header+getRequiredBytes(),reserveEnd+sizeof_Header);
        if (device) device->memoryMonitor(newSize-sizeof_Header-allocEnd,true);
        reserveEnd = allocEnd = newSize-sizeof_Header;
        if (next) next->shrink(device);
      }

      size_t getRequiredBytes() const {
        return min(size_t(cur),reserveEnd);
      }

      size_t getBlockAllocatedBytes() const {
	return min(max(allocEnd,size_t(cur)),reserveEnd);
      }

      size_t getAllocatedBytes() const {
	return getBlockAllocatedBytes() + (next ? next->getAllocatedBytes() : 0);
      }

      size_t getReservedBytes() const {
	return reserveEnd + (next ? next->getReservedBytes() : 0);
      }

      size_t getFreeBytes() const {
	return max(allocEnd,size_t(cur))-cur;
      }

      void print() const {
        std::cout << "[" << cur << ", " << max(size_t(cur),allocEnd) << ", " << reserveEnd << "] ";
        if (next) next->print();
      }

    public:
      std::atomic<size_t> cur;        //!< current location of the allocator
      std::atomic<size_t> allocEnd;   //!< end of the allocated memory region
      std::atomic<size_t> reserveEnd; //!< end of the reserved memory region
      Block* next;               //!< pointer to next block in list
      char align[maxAlignment-4*sizeof(size_t)]; //!< align data to maxAlignment
      char data[1];              //!< here starts memory to use for allocations
    };


  private:
    MemoryMonitorInterface* device;
    SpinLock mutex;
    size_t slotMask;
    std::atomic<Block*> threadUsedBlocks[MAX_THREAD_USED_BLOCK_SLOTS];
    std::atomic<Block*> usedBlocks;
    std::atomic<Block*> freeBlocks;

    std::atomic<Block*> threadBlocks[MAX_THREAD_USED_BLOCK_SLOTS];
    SpinLock slotMutex[MAX_THREAD_USED_BLOCK_SLOTS];

    size_t growSize;
    size_t bytesUsed;            //!< bumber of total bytes used

    ThreadLocalData<ThreadLocal> thread_local_allocators; //!< thread local allocators
    ThreadLocalData<ThreadLocal2> thread_local_allocators2; //!< thread local allocators
  };
}
