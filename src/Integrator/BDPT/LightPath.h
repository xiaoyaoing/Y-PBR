#pragma  once

#include "PathVertex.h"
#include "scene.hpp"
#include "Sampler/Sampler.hpp"

class LightPath {
    std::unique_ptr<PathVertex[]> _vertexs;
    int maxlength;
    bool adjoint;
    int _length;
public :
    inline int getLength() const {
        return _length;
    }

    const PathVertex &operator[](int i) const {
        return _vertexs.operator[](i);
    }

    void startCameraPath(const Camera *camera, ivec2 point);

    void startLightPath(const Light *light, Float lightPdf);

    void tracePath(const Scene &scene, Sampler &sampler, int traceMaxLength = -1);

    LightPath(int maxlength) : _length(0), maxlength(maxlength), adjoint(false),
                               _vertexs(new PathVertex[maxlength + 2]) {}

    static Spectrum
    connectCameraBDPT(const Scene &scene, const Camera *camera, Sampler &sampler, const LightPath &lightPath, int l,
                      ivec2 &pixel);

    static Spectrum
    connectLightBDPT(const Scene &scene, const Light *light, Sampler &sampler, const LightPath &cameraPath, int c,
                     Float lightPdf);

    static Spectrum connectBDPT(const Scene &scene, const LightPath &lightPath, int l, const LightPath &cameraPath,
                                int c);

    static inline Ray generateRay(const PathVertex &a, const PathVertex &b) {
        vec3 dir = b.pos() - a.pos();
        auto l = length(dir);
        dir /= l;
        return Ray(a.pos(), dir, 1e-4f, l - 1e-4f);
    }

};