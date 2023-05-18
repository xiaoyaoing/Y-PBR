#include "BdptTracer.h"

Spectrum BdptTracer::traceSample(ivec2 pixel, Sampler &sampler) {
    cameraPath->startCameraPath(camera, pixel);
    cameraPath->tracePath(scene, sampler);

    float lightPdf;
    auto lightIdx = lightDistrib->SampleDiscrete(sampler.getNext1D(), &lightPdf);
    const auto &light = scene.lights[lightIdx];
    auto lightSample = light->sampleDirect(sampler.getNext2D(), sampler.getNext2D());
    lightPath->startLightPath(light.get(), lightPdf);
    lightPath->tracePath(scene,sampler);

    int cameraLength = cameraPath->getLength();
    int lightLength = lightPath->getLength();

    Spectrum L(0);

    for (int l = 1; l <= lightLength; ++l){
        int upperBound = std::min(maxBounces -l +1, cameraLength);
        for (int c = 1; c <= upperBound; ++c) {
            if(l!=2 || c!=1)
                continue;
                if (!cameraPath->operator[](c-1).canConnect() || !lightPath->operator[](l-1).canConnect())
                    continue;
            if (c == 1) {
                ivec2 pRaster;
                Spectrum s = LightPath::connectCameraBDPT(scene, camera, sampler, *lightPath, l, pRaster);
                image->addPixel(pixel.x, pixel.y, s);
                if (imagePramid)
                    imagePramid->addPixel(l, c, pRaster, s);
            } else if (l == 1) {
                Spectrum s = LightPath::connectLightBDPT(scene, light.get(), sampler, *cameraPath, c,lightPdf);
                if(luminace(s)<0){
                    int k =1;
                }
                L += s;
                if (imagePramid)
                    imagePramid->addPixel(l, c, pixel, s);
            } else {
                Spectrum s = LightPath::connectBDPT(scene, *lightPath, l, *cameraPath, c);
                if(luminace(s)<0){
                    int k =1;
                }
                L += s;
            }
        }
    }
    return L;
}

BdptTracer::BdptTracer(const Scene &scene,const Distribution1D * distrib, const Camera *camera, ImagePramId *imagePramId,
                       int maxBounces) : scene(scene), camera(camera), image(camera->image.get()),
                                         lightDistrib(distrib),
                                         imagePramid(imagePramId),
                                         lightPath(std::make_unique<LightPath>(maxBounces + 1)),
                                         cameraPath(std::make_unique<LightPath>(maxBounces + 1)),
                                         maxBounces(maxBounces){

}



