#include "BDPTIntegrator.hpp"
#include "Sampler/LightDistrib.hpp"
#include "Bsdfs/Reflection.hpp"
void BDPTIntegrator::process(const Scene &scene, Sampler &sampler) {
    lightDistrib = CreateLightSampleDistribution(std::string("uniform"), scene);

}

vec3 BDPTIntegrator::integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const {

    return vec3();
}

void BDPTIntegrator::render(const Scene &scene) const {
    SamplerIntegrator::render(scene);
}


int generateLightPath(const Scene &  scene,const Distribution1D * lightDistrib, Sampler &sampler, int maxDepth, PathVertex * path) {
    float lightPdf;
    auto lightIdx =  lightDistrib->SampleDiscrete(sampler.getNext1D(),&lightPdf);
    const auto & light = scene.lights[lightIdx];
    auto lightSample =  light->sampleDirect(sampler.getNext2D(),sampler.getNext2D());
    path[0]  = PathVertex(light.get(),lightSample.radiance,lightSample.lightDirPdf * lightSample.lightPosPdf * lightPdf);
    return randomWalk(scene,lightSample.ray,sampler, maxDepth, path + 1);
}

int generateCameraPath(const Scene & scene,const Camera * camera,vec2 point, Sampler & sampler,int maxDepth,PathVertex * path){
    auto ray  = camera->sampleRay(point.x,point.y,sampler.getNext2D());
    path[0] = PathVertex(camera,point);
    return randomWalk(scene, ray,sampler, maxDepth, path + 1);
}

//only surface
int randomWalk(const Scene &scene, Ray &ray, Sampler &sampler, int maxDepth, PathVertex *path,bool adjoint,Float pdf) {
    auto its = scene.intersect(ray);
    int depth = 0;
    Float pdfFwd = pdf ,pdfBack;
    Spectrum  beta;
    while(its.has_value() && depth<maxDepth){
        PathVertex & vertex = path[depth],& prev = path[depth-1];

        SurfaceEvent event = makeLocalScatterEvent(& its.value());
        Spectrum  f = its->bsdf->sampleF(event,sampler.getNext2D(),adjoint);
        beta *= f/pdfFwd;
        if(russian(depth,sampler,beta))
            break;
        ray = event.sctterRay();

        path[depth++] = PathVertex(event,beta,pdfFwd,prev);
        pdfFwd = event.pdf;
        pdfBack = event.its->bsdf->Pdf(event.makeFlipQuery());
        prev.pdfBack = pdfBack;
    }
    return depth;
}

