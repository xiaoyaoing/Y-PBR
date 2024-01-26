#include "Integrator.hpp"
#include "Bsdfs/Reflection.hpp"
#include "Common/ProgressReporter.h"
#include "Sampler/Sampler.hpp"
#include "Common/Parallel.h"
#include "Camera/Camera.hpp"
#include "PathIntegrator.hpp"
#include <thread>
#include <iostream>

#include "Texture/BitMapTexture.hpp"

#include "Mediums/Medium.hpp"

static bool sampleBSDF  = true;
static bool sampleLgiht = true;

//Integrator::Integrator(Json j) {
//
//}

void SamplerIntegrator::render(const Scene& scene) {
    auto  tileSize     = scene.options.tileSize;
    ivec2 renderBounds = _camera->image->resoulation();
    int   width        = _camera->image->width();
    int   height       = _camera->image->height();
    ivec2 numTiles{(renderBounds.x + tileSize - 1) / tileSize, (renderBounds.y + tileSize - 1) / tileSize};

    int num_threads = std::thread::hardware_concurrency();
    parallel_init(num_threads);

    int spp     = scene.options.spp;
    int sppStep = scene.options.sppStep;

    ProgressReporter reporter(numTiles.x * numTiles.y);
    parallel_for([&](const vec2& tile) {
        int x0 = tile[0] * tileSize;
        int x1 = std::min(x0 + tileSize, width);
        int y0 = tile[1] * tileSize;
        int y1 = std::min(y0 + tileSize, height);

        int                      seed        = x0 * width + y0;
        std::unique_ptr<Sampler> tileSampler = _sampler->clone(seed);

        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                tileSampler->startPixel(ivec2(x, y));
                do {
                    {
                        Ray      ray      = _camera->sampleRay(x, y, tileSampler->getNext2D());
                        Spectrum radiance = integrate(ray, scene, *tileSampler);
                        _camera->image->addPixel(x, y, radiance, true);
                    }
                } while (tileSampler->startNextSample());
            }
        }
        reporter.update(1);
    },
                 numTiles);
    _camera->image->save(scene.options.outputFileName, 1.f, scene.options.overwrite);

    parallel_cleanup();
}

void SamplerIntegrator::renderPixel(int x, int y) const {
}