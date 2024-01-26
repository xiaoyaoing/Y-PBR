#include "BdptTracer.h"

Spectrum BdptTracer::traceSample(ivec2 pixel, Sampler& sampler) {
    cameraPath->startCameraPath(camera, pixel);
    cameraPath->tracePath(scene, sampler);

    float       lightPdf;
    auto        lightIdx    = lightDistrib.SampleDiscrete(sampler.getNext1D(), &lightPdf);
    const auto& light       = scene.lights[lightIdx];
    auto        lightSample = light->sampleDirect(sampler.getNext2D(), sampler.getNext2D());
    lightPath->startLightPath(light.get(), lightPdf);
    lightPath->tracePath(scene, sampler);

    int cameraLength = cameraPath->getLength();
    int lightLength  = lightPath->getLength();

    Spectrum L(0);
    for (int l = 0; l <= lightLength; ++l) {
        int upperBound = std::min(maxBounces - l + 1, cameraLength);
        for (int c = 1; c <= upperBound; ++c) {
            if (l + c < 2)
                continue;
            Spectrum pathContributionWithoutMis;
            ivec2    pRaster;

            if (l == 0) {
                pathContributionWithoutMis = LightPath::cameraDirectLight(scene, *cameraPath, c);
                //  continue;
            } else {
                if (!cameraPath->operator[](c - 1).canConnect() || !lightPath->operator[](l - 1).canConnect()) {
                    continue;
                }

                if (c == 1) {
                    pathContributionWithoutMis = LightPath::connectCameraBDPT(scene, sampler, *lightPath, *cameraPath, l, pRaster);
                    //  continue;
                } else if (l == 1) {
                    pathContributionWithoutMis = LightPath::connectLightBDPT(scene, sampler, *lightPath, *cameraPath, c, lightPdf);
                    // continue;
                    //   pathContributionWithoutMis = LightPath::connectBDPT(scene, *lightPath, l, *cameraPath, c);

                } else {
                    pathContributionWithoutMis = LightPath::connectBDPT(scene, *lightPath, l, *cameraPath, c);
                    // continue;
                }
            }
            if (isBlack(pathContributionWithoutMis))
                continue;
            pathContributionWithoutMis *= LightPath::misWeight(*lightPath, l, *cameraPath, c, lightDistrib, lightIdxs);
            // L = Spectrum(LightPath::misWeight(*lightPath, l, *cameraPath, c, lightDistrib, lightIdxs));
            // pathContributionWithoutMis = L;
            //   pathContributionWithoutMis = Spectrum(1/LightPath::misWeight(*lightPath, l, *cameraPath, c, lightDistrib, lightIdxs));
            if (c == 1) {
                {
                    image->addPixel(pRaster.x, pRaster.y, pathContributionWithoutMis, false);
                    if (imagePramid)
                        imagePramid->addPixel(l, c, pRaster, pathContributionWithoutMis, false);
                    continue;
                }
            }
            L += pathContributionWithoutMis;
            if (imagePramid)
                imagePramid->addPixel(l, c, pixel, pathContributionWithoutMis, false);
        }
    }
    return L;
}

BdptTracer::BdptTracer(const Scene& scene, const Distribution1D& distrib, const std::map<const Light*, size_t>& lightIdxs, const Camera* camera, ImagePrid* imagePramId, int maxBounces) : scene(scene), camera(camera), image(camera->image.get()),
                                                                                                                                                                                           lightDistrib(distrib),
                                                                                                                                                                                           lightIdxs(lightIdxs),
                                                                                                                                                                                           imagePramid(imagePramId),
                                                                                                                                                                                           lightPath(std::make_unique<LightPath>(maxBounces + 1, false)),
                                                                                                                                                                                           cameraPath(std::make_unique<LightPath>(maxBounces + 1, false)),
                                                                                                                                                                                           maxBounces(maxBounces) {
}