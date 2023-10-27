#include "LightPath.h"
#include "Bsdfs/Reflection.hpp"

void LightPath::startCameraPath(const Camera *camera, ivec2 point) {
    _length = 1;
    _vertexs[0].initCamera(camera, point);
    adjoint = false;
}

void LightPath::startLightPath(const Light *light, Float lightPdf) {
    _length = 1;
    _vertexs[0].initLight(light, lightPdf);
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
        PathVertex &cur = _vertexs[_length - 1];
        if (!cur.sampleNext(scene, adjoint, state, prev, _vertexs[_length]))
            break;
        _length++;
    }
    toAreaMeasure();
}

bool onlyMisWeight = false;
bool onlyLight = false;

Spectrum
LightPath::connectCameraBDPT(const Scene &scene, Sampler &sampler,const LightPath &lightPath,  const LightPath &cameraPath, int l,
                             ivec2 &pixel) {

    auto camera = cameraPath[0]._sampler.camera;
    const auto &lightVertex = lightPath[l - 1];
    if (!lightVertex.canConnect())
        return Spectrum(0);
    PositionAndDirectionSample sample;

    if (!camera->sampleLi(lightVertex.pos(), &pixel, sampler.getNext2D(), sample))
        return Spectrum(0, 0, 0);

    auto cameraVertex = PathVertex(camera, sample);

    auto ray = generateRay(cameraVertex, lightVertex);

    auto s=  lightVertex.beta * cameraVertex.beta * lightVertex.eval(cameraVertex, true);

    if(isBlack(s))
        return Spectrum(0);

    auto tr = evalShadowDirect(scene, ray, nullptr);
    s *= tr;

    if(isBlack(s))
        return Spectrum(0);

    Float weight = misWeight(lightPath,l,cameraPath,1);
    if(onlyMisWeight)
        return  Spectrum(weight);

    if(onlyLight)
        return s;
    return s * weight;
}

Spectrum
LightPath::connectBDPT(const Scene &scene, const LightPath &lightPath, int l, const LightPath &cameraPath, int c) {

    const PathVertex &lightVertex = lightPath[l - 1];
    const PathVertex &cameraVertex = cameraPath[c - 1];

    auto s = lightVertex.eval(cameraVertex, true) * cameraVertex.eval(lightVertex, false) * lightVertex.beta *
             cameraVertex.beta
             /distance2(lightVertex.pos(), cameraVertex.pos());

    if(isBlack(s))
        return s;

    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = evalShadowDirect(scene, ray, nullptr);
    if(isBlack(s * tr))
        return Spectrum(0);
    Float weight = misWeight(lightPath,l,cameraPath,c);
    if(onlyMisWeight)
        return  Spectrum(weight);
    if(onlyLight)
        return s * tr;
    return s * tr * weight;

}

Spectrum
LightPath::connectLightBDPT(const Scene &scene, Sampler &sampler,const LightPath &lightPath,  const LightPath &cameraPath, int c,
                            Float lightPdf) {
    auto light = lightPath[0]._sampler.light;
    const auto &cameraVertex = cameraPath[c - 1];
    if (!cameraVertex.canConnect())
        return Spectrum(0);
    auto sample = light->sampleLi(cameraVertex.pos(), sampler.getNext2D());
    if (isBlack(sample.weight) || sample.posPdf == 0)
        return Spectrum(0);

    auto lightVertex = PathVertex(light, sample, lightPdf);
    auto ray = generateRay(cameraVertex, lightVertex);
    auto tr = evalShadowDirect(scene, ray, nullptr);
    if (isBlack(tr))
        return Spectrum(0, 0, 0);
    auto res = lightVertex.beta * cameraVertex.beta * cameraVertex.eval(lightVertex, false);
    if(isBlack(res))
        return Spectrum(0);
    //return  Spectrum(misWeight(lightPath,1,cameraPath,c));
    Float weight = misWeight(lightPath,1,cameraPath,c);
    if(onlyMisWeight)
        return  Spectrum(weight);
    if(onlyLight)
        return res;
    return res * weight;
}

Float LightPath::misWeight(const LightPath &lightPath, int l, const LightPath &cameraPath, int c) {

//    return 1.f;
 //   return (1.f/(l+c));
    Float *pdfForward = reinterpret_cast<float *>(alloca((l + c) * sizeof(float)));
    Float *pdfBackward = reinterpret_cast<float *>(alloca((l + c) * sizeof(float)));
    Float *r = reinterpret_cast<float *>(alloca((l + c) * sizeof(float)));
    bool *connectable = reinterpret_cast<bool *>(alloca((l + c) * sizeof(bool)));

    //todo 判断相机路径末尾和光线路径末尾能否连接

    for (int i = 0; i < l; ++i) {
        pdfForward[i] = lightPath[i].pdfFwd;
        pdfBackward[i] = lightPath[i].pdfBack;
        connectable[i] = !lightPath[i].isDelta();
    }
    for (int i = 0; i < c; ++i) {
        pdfForward[l + c - (i + 1)] = cameraPath[i].pdfFwd;
        pdfBackward[l + c - (i + 1)] = cameraPath[i].pdfBack;
        connectable[l + c - (i + 1)] = !cameraPath[i].isDelta();
    }


    /// 跟其他可能构成同样长度的路径来比
    /// ps/sum(p0,p1,p2...pn)
    /// pi =  pf_0 * pf_1 * ... *  pf_i-1  * pb_i * ... pb_(n-1)
    /// mis weight : w_s
    /// let ri = pi / ps;
    /// w_s = 1/ (r0 + r1+r2+...+1+r_(s+1)+...r(n-1))
    /// ri = {
    ///         1  if i==s,
    ///         pb_i / pf(i) * r_i+1 if i <
    ///          pf_(i-1) / pb(i-1) * r_(i-1) if i <s
    /// }

    /// 相机路径和光子路径的最后一个节点的pdfBack没有计算
    const PathVertex *lightEnd = &lightPath[l - 1],
            *cameraEnd = &cameraPath[c - 1],
            *lightMinus = l > 1 ? &lightPath[l - 2] : nullptr,
            *cameraMinus = c > 1 ? &cameraPath[c - 2] : nullptr;

    // cameraEndPdf lightEndPdf cameraMinusPdf lightMinusPdf;
    pdfBackward[l] = lightEnd->pdf(lightMinus, *cameraEnd);
    pdfBackward[l-1] = cameraEnd->pdf(cameraMinus, *lightEnd);
    if (cameraMinus)
        pdfBackward[l + 1] = cameraEnd->pdf(lightEnd, *cameraMinus);
    if (lightMinus)
        pdfBackward[l - 2] = lightEnd->pdf(cameraEnd, *lightMinus);

    auto remap0 = [](Float f) -> Float { return f != 0 ? f : 1; };


    Float ri = 1, sumR = 1;
    for (int i = l ; i < l + c -1 ; i++) {


        ri *= remap0(pdfBackward[i]) / remap0(pdfForward[i]) ;

        if (connectable[i] && connectable[i + 1])
            sumR += ri;
    }
    ri=1;
    for (int i = l - 1; i >= 1; i--) {
        ri *= remap0(pdfBackward[i]) / remap0(pdfForward[i]);
        if (connectable[i] && connectable[i - 1])
            sumR += ri;
    }

    return 1.f / sumR;

}

void LightPath::toAreaMeasure() {
   // return ;
    for(int i=1;i<_length;i++){
        if(_vertexs[i-1].isDelta())
            continue;
        _vertexs[i].pdfFwd /= distance2(_vertexs[i].pos(),_vertexs[i-1].pos());
        auto d= distance2(_vertexs[i].pos(),_vertexs[i-1].pos());
        if(i==2 && _vertexs[i].pdfFwd <=0.003){
            int k  = 1;
        }
        if(_vertexs[i].isSurface())
            _vertexs[i].pdfFwd *= _vertexs[i-1].cosFactor(_vertexs[i]);
        //_vertexs[i].pdfFwd = 1/ distance2(_vertexs[i].pos(),_vertexs[i-1].pos());
        if(isinf(_vertexs[i].pdfFwd)){
            int k  = 1;
        }
    }
    for(int i= _length-3;i>=0;i--){
        if (_vertexs[i + 1].isDelta())
            continue;
        _vertexs[i].pdfBack /= distance2(_vertexs[i].pos(),_vertexs[i+1].pos());
        if (_vertexs[i].isSurface())
            _vertexs[i].pdfBack *= _vertexs[i].cosFactor(_vertexs[i+1]);
    }
}
