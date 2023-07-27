#include "LightPath.h"
#include "Bsdfs/Reflection.hpp"
void LightPath::startCameraPath(const Camera *camera, ivec2 point) {
    _length = 1;
    _vertexs[0].initCamera(camera,point);
    adjoint = false;
}

void LightPath::startLightPath(const Light *light, Float lightPdf) {
    _length = 1;
    _vertexs[0].initLight(light,lightPdf);
 //   _vertexs[0] = PathVertex(light, lightPdf);
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
        return Spectrum(0,0,0);
    //return Spectrum(1);
    auto cameraVertex = PathVertex(camera,sample);

    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = evalShadowDirect(scene, ray, nullptr);


    auto res = tr *   lightVertex.beta   * cameraVertex.beta * lightVertex.eval(cameraVertex, true);
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

Float LightPath::misWeight(const LightPath &lightPath, int l, const LightPath &cameraPath, int c) {
    Float *pdfForward           = reinterpret_cast<float *>(alloca((l+c)*sizeof(float)));
    Float *pdfBackward          = reinterpret_cast<float *>(alloca((l+c)*sizeof(float)));
    Float *r          = reinterpret_cast<float *>(alloca((l+c)*sizeof(float)));
    bool  *connectable          = reinterpret_cast<bool  *>(alloca((l+c)*sizeof(bool)));

    //todo 判断相机路径末尾和光线路径末尾能否连接

    for (int i = 0; i < l; ++i) {
        pdfForward [i] = lightPath[i].pdfFwd;
        pdfBackward[i] = lightPath[i].pdfBack;
        connectable[i] = lightPath[i].canConnect();
    }
    for (int i = 0; i < c; ++i) {
        pdfForward [l + c - (i + 1)] = cameraPath[i].pdfFwd;
        pdfBackward[l + c - (i + 1)] = cameraPath[i].pdfBack;
        connectable[l + c - (i + 1)] = cameraPath[i].canConnect();
    }


    /// 跟其他可能构成同样长度的路径来比
    /// ps/sum(p0,p1,p2...pn)
    /// pi =  pf_0 * pf_1 * ... *  pf_i-1  * pb_i * ... pb_(n-1)
    /// mis weight : w_s
    /// let ri = pi / ps;
    /// w_s = 1/ (r0 + r1+r2+...+1+r_(s+1)+...r(n-1))
    /// ri = {
    ///         1  if i==s,
    ///         pb_i / pf(i) * r_i+1 if i <s
    ///          pf_(i-1) / pb(i-1) * r_(i-1) if i <s
    /// }

    /// 相机路径和光子路径的最后一个节点的pdfBack没有计算
    

    Float ri = 1,sumR = 1;
    for(int i=l-1;i>=0;i--){
        ri *= pdfBackward[i] / pdfForward[i];
        sumR += ri;
    }
    ri = 1;
    for(int i=0;i<l;i++){
        ri *=   pdfForward[i-1] / pdfBackward[i-1];
        sumR += ri;
    }
    return 1.f/sumR;

}
