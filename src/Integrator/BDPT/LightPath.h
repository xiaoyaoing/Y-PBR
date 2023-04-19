#pragma  once
#include "PathVertex.h"
#include "scene.hpp"
#include "Sampler/Sampler.hpp"
class LightPath{
    std::unique_ptr<PathVertex[]> _vertexs;
    int maxlength;
    bool adjoint;
public :
    void startCameraPath(const Camera *camera, ivec2 point);
    void startLightPath(const Light *light, Float lightPdf);
    void tracePath(const Scene & scene,Sampler & sampler,int traceMaxLength = -1);
    LightPath(int maxlength) : length(0), maxlength(maxlength),adjoint(false),_vertexs(new PathVertex[maxlength+2]){}
    PathVertex & operator [](int idx)const{
        return _vertexs.operator[](idx);
    }
    static Spectrum
    connectCameraBDPT(const Scene &scene, const Camera *camera, Sampler &sampler, const LightPath &lightPath, int l,
                      ivec2 &pixel);
    static Spectrum
    connectLightBDPT(const Scene &scene, const Light *light, Sampler &sampler, const LightPath &cameraPath, int c,
                     );
    static Spectrum  connectBDPT(const Scene &scene,  const LightPath & lightPath, int l, const LightPath & cameraPath,
                       int c);
    static inline Ray  generateRay(const PathVertex & a,const PathVertex & b){
        auto dir = b.pos() -a.pos();
        auto l = dir.length();
        dir/=l;
        return Ray(a.pos(),dir,Constant::EPSILON,l-Constant::EPSILON);
    }
    int length;

};