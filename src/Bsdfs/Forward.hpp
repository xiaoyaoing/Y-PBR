#include "Reflection.hpp"

//Forward BSDF
//Light will ignore this BSDF
//Used for volume rendering
class ForwardBSDF : public BSDF {
public:
    ForwardBSDF() : BSDF(BSDF_FORWARD) {}
    Float Pdf(const SurfaceEvent& event) const override;

protected:
    Spectrum sampleF(SurfaceEvent& event, const vec2& u) const override;

    Spectrum f(const SurfaceEvent& event) const override;
};