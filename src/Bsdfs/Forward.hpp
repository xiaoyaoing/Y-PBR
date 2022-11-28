#include "Reflection.hpp"

class ForwardBSDF : public BSDF{
public:
    ForwardBSDF(): BSDF(BSDF_FORWARD) {}
    Float Pdf(const SurfaceEvent & event) const override;

    void LogInfo( ) const override;

protected:
    Spectrum sampleF(SurfaceEvent & event, const vec2 & u) const override;

    Spectrum f(const SurfaceEvent & event) const override;
};