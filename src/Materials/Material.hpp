#include "Bsdfs/Reflection.hpp"
#include "Ray/Intersection.hpp"
#include "Common/Texture.hpp"
class Material{
public:
    // Material Interface
    virtual void ComputeScatteringFunctions(SurfaceIntersection *si);
//                                            MemoryArena &arena,
//                                            TransportMode mode,
//                                            bool allowMultipleLobes)
//

    virtual ~Material();
    static void Bump(const std::shared_ptr<Texture<Float>> &d,SurfaceIntersection *si);
};

