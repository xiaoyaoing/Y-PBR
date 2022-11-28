#include "Forward.hpp"

Float ForwardBSDF::Pdf(const SurfaceEvent & event) const {
    return 0;
}

void ForwardBSDF::LogInfo( ) const {

}

Spectrum ForwardBSDF::sampleF(SurfaceEvent & event, const vec2 & u) const {
    return Spectrum();
}

Spectrum ForwardBSDF::f(const SurfaceEvent & event) const {
    return Spectrum();
}
