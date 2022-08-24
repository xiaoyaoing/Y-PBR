#include "Primitive.hpp"


Intersection Primitive::Sample(const Intersection & ref, const vec2 & u, Float * pdf) const {
    auto intr =  Sample(u,pdf);

    vec3 wi = intr.p - ref.p;
    if (length2(wi) == 0)
        *pdf = 0;
    else {
        wi = normalize(wi);
        // Convert from area measure, as returned by the Sample() call
        // above, to solid angle measure.
        *pdf *= length2(ref.p-intr.p) / abs(dot(intr.getNormal(), -wi));
        if (std::isinf(*pdf)) *pdf = 0.f;
    }
    return intr;
}


