#include "BoundingBox.hpp"

bool Bounds3::IntersectP(const Ray& ray, const vec3& invDir, const int* dirIsNeg) const {
    const Bounds3& bounds = *this;
    // Check for ray intersection against $x$ and $y$ slabs
    Float tMin  = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
    Float tMax  = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
    Float tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
    Float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

    // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
    tMax *= 1 + 2 * 1e-3f;
    tyMax *= 1 + 2 * 1e-3f;
    if (tMin > tyMax || tyMin > tMax) return false;
    if (tyMin > tMin) tMin = tyMin;
    if (tyMax < tMax) tMax = tyMax;

    // Check for ray intersection against $z$ slab
    Float tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
    Float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

    // Update _tzMax_ to ensure robust bounds intersection
    tzMax *= 1 + 2 * 1e-3f;
    if (tMin > tzMax || tzMin > tMax) return false;
    if (tzMin > tMin) tMin = tzMin;
    if (tzMax < tMax) tMax = tzMax;
    return (tMin < ray.farT) && (tMax > 0);
}