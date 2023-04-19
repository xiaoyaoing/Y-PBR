#include "BDPTIntegrator.hpp"
#include "Sampler/LightDistrib.hpp"
#include "Bsdfs/Reflection.hpp"
#include "Common/Parallel.h"
#include "Common/ProgressReporter.h"
#include <thread>
void BDPTIntegrator::process(const Scene &scene, Sampler &sampler) {
    lightDistrib = CreateLightSampleDistribution(std::string("uniform"), scene);
}

vec3 BDPTIntegrator::integrate(const Ray &ray, const Scene &scene, Sampler &sampler) const {

    return vec3();
}

void BDPTIntegrator::render(const Scene &scene)  {
    auto tileSize = scene.options.tileSize;
    ivec2 renderBounds = _camera->image->resoulation();
    int width = _camera->image->width();
    int height = _camera->image->height();
    ivec2 numTiles{(renderBounds.x + tileSize - 1) / tileSize, (renderBounds.y + tileSize - 1) / tileSize};
    for(int i=0 ;i<numTiles.x * numTiles.y;i++)
    {
        _tracers.emplace_back(BdptTracer(scene,_camera.get(),imagePramid,scene.options.maxBounces));
    }
    int num_threads = std::thread::hardware_concurrency();
    parallel_init(num_threads);

    int spp = scene.options.spp;
    int sppStep = scene.options.sppStep;

    ProgressReporter reporter(numTiles.x * numTiles.y);
    parallel_for([&](const vec2 &tile) {

        int x0 = tile[0] * tileSize;
        int x1 = std::min(x0 + tileSize, width);
        int y0 = tile[1] * tileSize;
        int y1 = std::min(y0 + tileSize, height);

        std::unique_ptr<Sampler> tileSampler = _sampler->clone();
        tileSampler->setSeed(tile.y * renderBounds.x + tile.x);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                for (int s = 0; s < spp; s++) {
                    {
                    }
                }
            }
        }
        reporter.update(1);
    }, numTiles);
    saveOutPuts(scene.options.outputFileName);
    parallel_cleanup();
}

void BDPTIntegrator::saveOutPuts(const std::string & fileName) {
    _camera->image->postProgress();
    _camera->image->savePNG(fileName);
    if(imagePramid){
       imagePramid->saveOutPut(fileName);
    }
}


//int generateLightPath(const Scene &  scene,const Distribution1D * lightDistrib, Sampler &sampler, int maxDepth, PathVertex * path) {
//    float lightPdf;
//    auto lightIdx =  lightDistrib->SampleDiscrete(sampler.getNext1D(),&lightPdf);
//    const auto & light = scene.lights[lightIdx];
//    auto lightSample =  light->sampleDirect(sampler.getNext2D(),sampler.getNext2D());
//    path[0]  = PathVertex(light.get(), lightSample.radiance, lightSample.dirPdf * lightSample.posPdf * lightPdf);
//    return randomWalk(scene,lightSample.ray,sampler, maxDepth, path + 1,true,1);
//}
//
//int generateCameraPath(const Scene & scene,const Camera * camera,vec2 point, Sampler & sampler,int maxDepth,PathVertex * path){
//    auto ray  = camera->sampleRay(point.x,point.y,sampler.getNext2D());
//    path[0] = PathVertex(camera,point);
//    return randomWalk(scene, ray,sampler, maxDepth, path + 1,false,1);
//}

////only surface
//int randomWalk(const Scene &scene, Ray &ray, Sampler &sampler, int maxDepth, PathVertex *path,bool adjoint,Float pdf) {
//    auto its = scene.intersect(ray);
//    int depth = 0;
//    Float pdfFwd = pdf ,pdfBack;
//    Spectrum  beta;
//    while(its.has_value() && depth<maxDepth){
//        PathVertex & vertex = path[depth],& prev = path[depth-1];
//
//        SurfaceEvent event = makeLocalScatterEvent(& its.value());
//        Spectrum  f = its->bsdf->sampleF(event,sampler.getNext2D(),adjoint);
//        beta *= f/pdfFwd;
//        if(russian(depth,sampler,beta))
//            break;
//        ray = event.sctterRay();
//
//        path[depth++] = PathVertex(event,beta,pdfFwd,prev);
//        pdfFwd = event.pdf;
//        pdfBack = event.its->bsdf->Pdf(event.makeFlipQuery());
//        prev.pdfBack = pdfBack;
//    }
//    return depth;
//}

void connectPath(PathVertex *lightPath, int ln, PathVertex *cameraPath, int cn, float *misWeight) {

}

//Spectrum
//connectPath(const Scene &scene, const Camera *camera, const PathVertex *lightPath, int ln, const PathVertex *cameraPath,
//            int cn, float *misWeight, ivec2 *pRaster) {
//    Spectrum  L;
//    PathVertex sampled;
//    if(ln == 0){
//        const PathVertex & pt = cameraPath[cn-1];
//        if(pt.isLight()){
//            //L = pt.
//        }
//    }
//    else if(cn == 1){
//        const PathVertex & lightEnd = lightPath[ln-1];
//        if(lightEnd.canConnect()){
//            const PathVertex & cameraVertex = cameraPath[0];
//            Float pdf;
//            Spectrum cameraS = camera->evalDirection(lightEnd.pos(), pRaster, &pdf);
//            if(isBlack(cameraS)){
//                return Spectrum(0);
//            }
//            Spectrum s = lightEnd.beta * lightEnd.eval(cameraVertex,true) * cameraVertex.beta * cameraS / pdf;
//            *misWeight = 1.f;
//            return s* *misWeight;
//        }
//    }
//}




