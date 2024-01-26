#pragma once

#include "math.hpp"

#include <mutex>
#include <functional>
#include <atomic>

// From https://github.com/mmp/pbrt-v3/blob/master/src/core/parallel.h
extern thread_local int ThreadIndex;

void parallel_for(const std::function<void(int64_t)>& func, int64_t count, int64_t chunk_size = 1);
void parallel_for(std::function<void(ivec2)> func, const ivec2 count);

void parallel_init(int num_threads);
void parallel_cleanup();