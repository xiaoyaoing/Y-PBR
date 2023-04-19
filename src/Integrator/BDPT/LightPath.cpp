#include "LightPath.h"

void LightPath::startCameraPath(const Camera *camera, ivec2 point) {
    length = 0;
    _vertexs[0] = PathVertex(camera, point);
    adjoint = false;
}

void LightPath::startLightPath(const Light *light, Float lightPdf) {
    length = 1;
    _vertexs[0] = PathVertex(light, lightPdf);
    adjoint = true;
}

void LightPath::tracePath(const Scene &scene, Sampler &sampler, int traceMaxLength) {
    PathState state(sampler);
    if (!_vertexs[0].sampleRootVertex(state))
        return;
    if (traceMaxLength == -1) traceMaxLength = maxlength;
    else traceMaxLength = std::min(traceMaxLength, maxlength);
    while (length < traceMaxLength) {
        PathVertex *prev = (length == 1 ? nullptr : &_vertexs[length - 1]);
        if (!_vertexs[length].sampleNext(scene, adjoint, state, prev, _vertexs[length]))
            break;
        length++;
    }
}

Spectrum
LightPath::connectCameraBDPT(const Scene &scene, const Camera *camera, Sampler &sampler, const LightPath &lightPath,
                             int l,
                             ivec2 &pixel) {
    const auto &lightVertex = lightPath[l - 1];
    if(!lightVertex.canConnect())
        return Spectrum();
    auto cameraVertex = PathVertex(camera, camera->sampleLi(lightVertex.pos(), &pixel, sampler.getNext2D()));
    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = TraceHelper::evalShadowDirect(scene, ray, nullptr);
    if (isBlack(tr))
        return Spectrum();
    return tr * cameraVertex.beta * lightVertex.beta * lightVertex.eval(cameraVertex, true);
}

Spectrum
LightPath::connectBDPT(const Scene &scene, const LightPath &lightPath, int l, const LightPath &cameraPath, int c) {
    const PathVertex &lightVertex = lightPath[l - 1];
    const PathVertex &cameraVertex = cameraPath[c - 1];
    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = TraceHelper::evalShadowDirect(scene, ray, nullptr);
    if (isBlack(tr))
        return Spectrum();
    if (lightVertex.canConnect() && cameraVertex.canConnect()) {
        return tr * lightVertex.eval(cameraVertex, true) * cameraVertex.eval(lightVertex, false) /
               distance2(lightVertex.pos(), cameraVertex.pos());
    }
}

Spectrum
LightPath::connectLightBDPT(const Scene &scene, const Light *light, Sampler &sampler, const LightPath &cameraPath,
                            int c) {
    const auto &cameraVertex = cameraPath[c - 1];
    if(!cameraVertex.canConnect())
        return Spectrum();
    auto lightVertex = PathVertex(light, light->sampleDirect(cameraVertex.pos(), sampler.getNext2D()));
    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = TraceHelper::evalShadowDirect(scene, ray, nullptr);
    if (isBlack(tr))
        return Spectrum();
    return tr * cameraVertex.beta * lightVertex.beta * cameraVertex.eval(lightVertex, false);
}
