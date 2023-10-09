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
//           if(l!=3)
//                continue;
            if (!cameraPath->operator[](c-1).canConnect() || !lightPath->operator[](l-1).canConnect())
            {
              //  imagePramid->addPixel(l, c, pixel, Spectrum(1,0,0));
                continue;
            }
            if (c == 1) {
                ivec2 pRaster;
                Spectrum s = LightPath::connectCameraBDPT(scene, sampler, *lightPath,*cameraPath,l, pRaster);
                image->addPixel(pRaster.x,pRaster.y, s);
                if (imagePramid)
                    imagePramid->addPixel(l, c, pRaster, s);
            } else if (l == 1) {
                Spectrum s = LightPath::connectLightBDPT(scene,  sampler,*lightPath,*cameraPath,c,lightPdf);
                L += s;
                if (imagePramid)
                    imagePramid->addPixel(l, c, pixel, s);
            } else {
                Spectrum s = LightPath::connectBDPT(scene, *lightPath, l, *cameraPath, c);
              //  s = vec3(1);
                if (imagePramid)
                    imagePramid->addPixel(l, c, pixel,s);
                L += s;
            }
        }
    }
    return L;
}

BdptTracer::BdptTracer(const Scene &scene, const Distribution1D * distrib, const Camera *camera, ImagePrid *imagePramId,
                       int maxBounces) : scene(scene), camera(camera), image(camera->image.get()),
                                         lightDistrib(distrib),
                                         imagePramid(imagePramId),
                                         lightPath(std::make_unique<LightPath>(maxBounces + 1)),
                                         cameraPath(std::make_unique<LightPath>(maxBounces + 1)),
                                         maxBounces(maxBounces){

}



