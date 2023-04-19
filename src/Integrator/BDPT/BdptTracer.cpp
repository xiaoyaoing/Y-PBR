#include "BdptTracer.h"

Spectrum BdptTracer::traceSample(ivec2 pixel, Sampler &sampler) {
    cameraPath->startCameraPath(camera,pixel);
    cameraPath->tracePath(scene,sampler);

    float lightPdf;
    auto lightIdx = lightDistrib->SampleDiscrete(sampler.getNext1D(),&lightPdf);
    const auto & light = scene.lights[lightIdx];
    auto lightSample =  light->sampleDirect(sampler.getNext2D(),sampler.getNext2D());
    lightPath->startLightPath(light.get(),lightPdf);

    int cameraLength = cameraPath->length;
    int lightLength = cameraPath->length;

    Spectrum L;
    for(int l = 0;l<lightLength;++l)
    {
        for (int c= 0; c < cameraLength; ++c) {
            if(c==1){
                ivec2 pRaster;
                Spectrum s = LightPath::connectCameraBDPT(scene,camera,sampler,*lightPath,l,pRaster);
                image->addPixel(pixel.x,pixel.y,s);
                if(imagePramid)
                    imagePramid->addPixel(l,c,pRaster,s);
            }
            else if(l==1) {
                Spectrum  s = LightPath::connectLightBDPT(scene,light.get(),sampler,*cameraPath,c);
                if(imagePramid)
                    imagePramid->addPixel(l,c,pixel,s);
            }
            else {
                Spectrum  s = LightPath::connectBDPT(scene,*lightPath,l,*cameraPath,c);
                L +=s;
                if(imagePramid)
                    imagePramid->addPixel(l,c,pixel,s);
            }
        }
    }
    return L;
}

BdptTracer::BdptTracer(const Scene &scene, const Camera *camera,  ImagePramId *imagePramId,
                       int maxBounces):scene(scene),camera(camera),image(camera->image.get()),imagePramid(imagePramId),
                                       lightPath(std::make_unique<LightPath>(maxBounces+1)),cameraPath(std::make_unique<LightPath>(maxBounces+1)){

}



