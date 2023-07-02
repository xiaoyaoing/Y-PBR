#include "LightPath.h"
#include "Bsdfs/Reflection.hpp"
void LightPath::startCameraPath(const Camera *camera, ivec2 point) {
    _length = 1;
    _vertexs[0] = PathVertex(camera, point);
    adjoint = false;
}

void LightPath::startLightPath(const Light *light, Float lightPdf) {
    _length = 1;
    _vertexs[0] = PathVertex(light, lightPdf);
    adjoint = true;
}

void LightPath::tracePath(const Scene &scene, Sampler &sampler, int traceMaxLength) {
    PathState state(sampler);
    if (!_vertexs[0].sampleRootVertex(state))
        return;
    if (traceMaxLength == -1) traceMaxLength = maxlength;
    else traceMaxLength = std::min(traceMaxLength, maxlength);
    while (_length < traceMaxLength) {
        PathVertex *prev = (_length == 1 ? nullptr : &_vertexs[_length - 2]);
        PathVertex & cur = _vertexs[_length-1];
        if (!cur.sampleNext(scene, adjoint, state, prev, _vertexs[_length]))
            break;
        _length++;
    }
}

Spectrum
LightPath::connectCameraBDPT(const Scene &scene, const Camera *camera, Sampler &sampler, const LightPath &lightPath,
                             int l,
                             ivec2 &pixel) {
    const auto &lightVertex = lightPath[l - 1];
    if(!lightVertex.canConnect())
        return Spectrum(0);
    PositionAndDirectionSample sample;
    if(!camera->sampleLi(lightVertex.pos(), &pixel, sampler.getNext2D(),sample))
        return Spectrum(1,0,0);

    auto cameraVertex = PathVertex(camera,sample);
    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = evalShadowDirect(scene, ray, nullptr);
   // if (isBlack(tr))
       // return Spectrum(0);
   // return lightVertex.eval(cameraVertex, true);
   // return lightVertex.beta;
   return lightVertex.eval(cameraVertex, true);
    auto res = tr *   lightVertex.beta  *  lightVertex.eval(cameraVertex, true) * cameraVertex.beta;
    return res;
    return lightVertex.eval(cameraVertex, true);
  //  return abs(lightVertex.pos())/10.f;
    return res;
    return cameraVertex.beta ;//* 10.f;
    if(l == 2 && distance2(lightVertex.pos(),lightPath[0].pos())<0.25 ){
        auto d = distance2(lightVertex.pos(),lightPath[0].pos());
        auto f= lightVertex.eval(cameraVertex, true);
        int k =1;
    }
    return  res;
}

Spectrum
LightPath::connectBDPT(const Scene &scene, const LightPath &lightPath, int l, const LightPath &cameraPath, int c) {
    const PathVertex &lightVertex = lightPath[l - 1];
    const PathVertex &cameraVertex = cameraPath[c - 1];
   // return lightVertex.beta * lightVertex.eval(cameraVertex, true) /  distance2(lightVertex.pos(), cameraVertex.pos())  ;
    auto s = lightVertex.eval(cameraVertex, true) * cameraVertex.eval(lightVertex, false) * lightVertex.beta * cameraVertex.beta /
             distance2(lightVertex.pos(), cameraVertex.pos());
    if(isBlack(s))
        return s;
    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = evalShadowDirect(scene, ray, nullptr);
    return s*tr;

}

Spectrum
LightPath::connectLightBDPT(const Scene &scene, const Light *light, Sampler &sampler, const LightPath &cameraPath,
                            int c,Float lightPdf) {
    const auto &cameraVertex = cameraPath[c - 1];
    if(!cameraVertex.canConnect())
        return Spectrum(0);
    auto sample  = light->sampleLi(cameraVertex.pos(),sampler.getNext2D());
    if(isBlack(sample.weight) || sample.posPdf ==0)
        return Spectrum(0);

    auto lightVertex = PathVertex(light, sample,lightPdf);
    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = evalShadowDirect(scene, ray, nullptr);
    //tr = Spectrum(1);
    if (isBlack(tr))
        return Spectrum(0,0,0);
    //return Spectrum(1);
    if(isBlack(lightVertex.beta * cameraVertex.beta) || luminace(lightVertex.beta * cameraVertex.beta)<1e-3f){
        int k = 1;
    }
    if(cameraVertex._sampler.bsdf->HasFlag(BXDFType(BSDF_SPECULAR | BSDF_REFLECTION))){
     //   return Spectrum(0);
    }
    auto res =  lightVertex.beta * cameraVertex.beta * cameraVertex.eval(lightVertex, false) ;
    if(hasNeg(res)){
        int k = 1;
    }
    return res;
}
