#include "PhotonMapper.hpp"
#include "Common/Parallel.h"
#include "Common/ProgressReporter.h"
#include "Bsdfs/Reflection.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "Sampler/LowDiscrepancy.hpp"
#include "Sampler/UniformSampler.h"
#include "SampleRecords/PositionAndDirectionSample.h"
#include "IO/ImageIO.hpp"
#include "TraceHelper.h"
#include <thread>
//sppm mainly learned from pbrt

struct SPPMPIxel {
    Float              radius     = 0;
    Float              lastRadius = 0;
    Spectrum           Ld         = Spectrum(0);
    Spectrum           flux       = Spectrum(0);
    std::atomic<Float> phi[3];
    vec2               pos;
    std::atomic<int>   curPhotonCount;
    Float              lastPhotonCount;
    float              rate;
    struct VisiblePoint {
        VisiblePoint() {}

        //        VisiblePoint(const vec3 & p, const vec3 & wo, const BSDF * bsdf, const Spectrum & beta)
        //                : p(p), wo(wo), bsdf(bsdf), beta(beta) {}
        SurfaceEvent* event = nullptr;
        Spectrum      beta  = Spectrum(0);

        ~VisiblePoint() {
            //delete event;
        }
    } vp;
};

struct SPPMListNode {
    SPPMPIxel*    pixel{nullptr};
    SPPMListNode* next{nullptr};
    bool          deleted{false};
    int           length() {
        auto node = this;
        int  ans  = 0;
        while (node && node->pixel) {
            ans++;
            node = node->next;
        }
        return ans;
    }

    ~SPPMListNode() {
        if (next) delete next;
    }
};

struct AtmoicVec3 {
    std::atomic<vec3> value;

    AtmoicVec3() {
        std::atomic_init(&value, vec3(0));
    }

    void add(const vec3& v) {
        value = std::atomic_load(&value) + v;
    }

    void divide(const Float a) {
        value = atomic_load(&value) / a;
    }

    void reset() {
        std::atomic_init(&value, vec3(0));
    }
};

static vec3 toGridPos(const Bounds3& bounds, const vec3& worldPos, const vec3& gridRes) {
    vec3 clampWorldPos = clamp(worldPos, bounds.pMin, bounds.pMax);
    auto t             = (clampWorldPos - bounds.pMin);
    auto b             = bounds.Diagonal();
    auto res           = t / b * gridRes;
    return res;
}

//convert 3d grid coordinate to linear index
static inline int hash(const ivec3& p, int hashSize) {
    return static_cast<unsigned int>((p.x * 73856093) ^ (p.y * 19349663) ^
                                     (p.z * 83492791)) %
           hashSize;
}

void PhotonMapper::render(const Scene& scene) {
    ivec2                        pixelBounds = _camera->image->resoulation();
    int                          pixelNum    = pixelBounds.x * pixelBounds.y;
    std::unique_ptr<SPPMPIxel[]> pixels(new SPPMPIxel[pixelNum]);
    for (int i = 0; i < pixelNum; i++) {
        pixels[i].radius     = initRadius;
        pixels[i].lastRadius = initRadius;
    }
    UniformSampler sampler;
    constexpr int  tileSize = 16;
    ivec2          nTiles((pixelBounds.x + tileSize - 1) / tileSize,
                 (pixelBounds.y + tileSize - 1) / tileSize);

    int num_threads = std::thread::hardware_concurrency();
    parallel_init(num_threads);
    //generate visible points
    ProgressReporter reporter(iterations);
    for (int iteration = 0; iteration < iterations; iteration++) {

        LOGI("Iteration {0} begin", iteration);

        //Generate pixel information
        parallel_for([&](ivec2 tileIndex) {
            std::unique_ptr<Sampler> tileSampler = sampler.clone((tileIndex.y * nTiles.x + tileIndex.x) * (iteration + 1));

            int x0 = tileIndex.x * tileSize;
            int x1 = std::min(x0 + tileSize - 1, pixelBounds.x - 1);
            int y0 = tileIndex.y * tileSize;
            int y1 = std::min(y0 + tileSize - 1, pixelBounds.y - 1);

            for (int y = y0; y <= y1; y++)
                for (int x = x0; x <= x1; x++) {

                    int   pixelIndex = y * pixelBounds.x + x;
                    auto& pixel      = pixels[pixelIndex];
                    pixel.pos        = {x, y};
                    Ray ray          = _camera->sampleRay(x, y, tileSampler->getNext2D());

                    bool     specularBounce = true;
                    Spectrum beta(1);
                    bool     isDielectric = false;

                    for (int bounce = 0; bounce < maxBounces; bounce++) {
                        std::optional<Intersection> its = scene.intersect(ray);
                        if (!its) {
                            for (auto light : scene.lights)
                                pixel.Ld += beta * light->Le(ray);
                            break;
                        }
                        //   pixel.Ld +=its->Ng;
                        // break;

                        if (bounce == 0 && its->bsdf->HasFlag(BSDF_SPECULAR)) {
                            isDielectric = true;
                        }
                        SurfaceEvent event = makeLocalScatterEvent(&its.value());
                        const BSDF*  bsdf  = its->bsdf;
                        if (specularBounce)
                            pixel.Ld += beta * its->Le(-ray.d);
                        // if(bounce!=0 && bsdf->HasFlag(BSDF_DIFFUSE))
                        pixel.Ld += beta * uniformSampleOneLight(event, scene, *tileSampler);
                        bool isDiffuse = bsdf->HasFlag(BSDF_DIFFUSE);
                        bool isGlossy  = bsdf->HasFlag(BSDF_GLOSSY);
                        //find visible {point}
                        if (isDiffuse || (isGlossy && bounce == maxBounces - 1)) {

                            pixel.vp.event      = new SurfaceEvent(event);
                            pixel.vp.event->its = new Intersection(*event.its);
                            pixel.vp.beta       = beta;
                            break;
                        }
                        if (bounce < maxBounces - 1) {
                            event.requestType = BSDF_ALL;
                            Spectrum f        = bsdf->sampleF(event, tileSampler->getNext2D(), false);
                            if (event.pdf == 0 || isBlack(f)) break;
                            specularBounce = (event.sampleType & BSDF_SPECULAR) != 0;
                            beta *= f / event.pdf;
                            if (luminace(beta) < 0.25) {
                                Float russinanProb = luminace(beta);
                                if (tileSampler->getNext1D() < russinanProb) {
                                    beta /= luminace(beta);
                                } else {
                                    break;
                                }
                            }
                            ray = event.sctterRay(event.toWorld(event.wi));
                            //   pixel.Ld += beta;
                            // break;
                        }
                    }
                    if (isDielectric) {
                    }
                }
        },
                     nTiles);

        //grid the visible points
        Bounds3 gridBounds;
        Float   maxRadius = 0;
        for (int i = 0; i < pixelNum; i++) {
            SPPMPIxel& pixel = pixels[i];
            if (isBlack(pixel.vp.beta))
                continue;
            Bounds3 vpBpunds(pixel.vp.event->its->p);
            vpBpunds   = Expand(vpBpunds, pixel.radius);
            gridBounds = Union(gridBounds, vpBpunds);
            maxRadius  = std::max(maxRadius, pixel.radius);
        }

        ivec3 gridRes;
        vec3  diag    = gridBounds.Diagonal();
        Float maxDiag = maxElement(diag);
        int   baseRes = static_cast<int>(maxDiag / maxRadius);
        for (int i = 0; i < 3; i++)
            gridRes[i] = std::max(static_cast<int>(baseRes * (diag[i] / maxDiag)), 1);

        std::vector<std::atomic<SPPMListNode*>> grids(pixelNum);
        parallel_for([&](int pixelIndex) {
            SPPMPIxel& pixel = pixels[pixelIndex];
            if (isBlack(pixel.vp.beta))
                return;
            vec3  pMin    = pixel.vp.event->its->p - vec3(pixel.radius);
            vec3  pMax    = pixel.vp.event->its->p + vec3(pixel.radius);
            ivec3 gridMin = toGridPos(gridBounds, pMin, gridRes);
            ivec3 gridMax = toGridPos(gridBounds, pMax, gridRes);
            for (int x = gridMin.x; x <= gridMax.x; x++)
                for (int y = gridMin.y; y <= gridMax.y; y++)
                    for (int z = gridMin.z; z <= gridMax.z; z++) {
                        int  gridIndex = hash(ivec3(x, y, z), pixelNum);
                        auto node      = new SPPMListNode();
                        node->pixel    = &pixel;
                        node->next     = grids[gridIndex].load();
                        while (!grids[gridIndex].compare_exchange_weak(
                            node->next, node))
                            ;
                    }
        },
                     pixelNum,
                     4096);
        LOGI("Visiable point generated");

        int              count = 0;
        std::vector<int> counts;
        for (auto& node : grids) {
            count += atomic_load(&node)->length();
            counts.push_back(atomic_load(&node)->length());
        }

        //trace Photons and accumlate contribution
        std::unique_ptr<Distribution1D> lightPowerDistrib = computeLightPowerDistrib(scene);
        parallel_for([&](int photonIndex) {
            uint64_t haltonIndex =
                static_cast<uint64_t>(iteration) * static_cast<uint64_t>(photonsPerIteration) +
                photonIndex + 1;
            int haltonDim = 0;

            bool drawLine = sampler.getNext1D() < 0.00000001;
            drawLine      = false;

            Spectrum lineColor = Spectrum(RadicalInverse(haltonDim, haltonIndex),
                                          RadicalInverse(haltonDim + 1, haltonIndex),
                                          RadicalInverse(haltonDim + 2, haltonIndex));
            haltonDim += 3;

            Float lightPdf;
            int   lightIndex = lightPowerDistrib->SampleDiscrete(RadicalInverse(haltonDim++, haltonIndex), &lightPdf);
            auto  light      = scene.lights[lightIndex];

            vec2 posSample(RadicalInverse(haltonDim, haltonIndex), RadicalInverse(haltonDim + 1, haltonIndex));
            vec2 dirSample(RadicalInverse(haltonDim + 2, haltonIndex), RadicalInverse(haltonDim + 3, haltonIndex));
            haltonDim += 5;
            PositionAndDirectionSample lightSample = light->sampleDirect(posSample, dirSample);

            if (lightSample.dirPdf == 0 || lightSample.posPdf == 0 || isBlack(lightSample.weight))
                return;
            Ray& photonRay = lightSample.ray;

            Spectrum beta = lightSample.weight * absDot(photonRay.d, lightSample.n) / (lightSample.dirPdf * lightSample.posPdf);
            //   beta = Spectrum(20.f);
            for (int bounce = 0; bounce < 8; bounce++) {
                std::optional<Intersection> its = scene.intersect(photonRay);
                if (!its.has_value()) {
                    break;
                }
                //Direct lighting has been considered, so contribution is calculated only when depth is greater than 1
                if (bounce > 0) {
                    vec3 photonP = its->p;
                    if (gridBounds.Contains(photonP)) {
                        vec3 gridPos   = toGridPos(gridBounds, photonP, gridRes);
                        int  gridIndex = hash(gridPos, pixelNum);
                        for (SPPMListNode* node = atomic_load(&grids[gridIndex]);
                             node != nullptr;
                             node = node->next) {

                            SPPMPIxel*               pixel = node->pixel;
                            SPPMPIxel::VisiblePoint& vp    = pixel->vp;
                            SurfaceEvent*            event = vp.event;
                            const vec3&              p     = event->its->p;
                            if (length2(p - photonP) > pixel->radius * pixel->radius)
                                continue;
                            event->wi        = event->toLocal(-photonRay.d);
                            Spectrum contrib = event->its->bsdf->f(*event, false) * beta;
                            for (int i = 0; i < 3; i++)
                                pixel->phi[i] = std::atomic_load(&pixel->phi[i]) + contrib[i];
                            ++pixel->curPhotonCount;
                        }
                    } else {
                    }
                }
                const BSDF*  photonBsdf = its->bsdf;
                SurfaceEvent event      = makeLocalScatterEvent(&its.value());
                event.requestType       = BSDF_ALL;
                Spectrum f              = photonBsdf->sampleF(event, vec2(RadicalInverse(haltonDim + 0, haltonIndex), RadicalInverse(haltonDim + 1, haltonIndex)), true);

                haltonDim += 2;
                Spectrum betaNew = beta * f / event.pdf;

                Float russProb = std::max(static_cast<Float>(0), 1 - luminace(betaNew) / luminace(beta));
                russProb       = 0;
                if (RadicalInverse(haltonDim++, haltonIndex) < russProb) {

                    break;
                }
                beta      = betaNew / (1 - russProb);
                vec3 wi   = event.toWorld(event.wi);
                photonRay = event.sctterRay(wi);
            }
        },
                     photonsPerIteration,
                     8192);

        for (auto& atomic_node : grids) {
            delete atomic_node.load();
        }

        LOGI("Trace photons and accumlate contribution completed");

        //update pixel
        std::atomic<float> rate        = 0;
        std::atomic<int>   updateCount = 0;
        std::atomic<float> curPhoton   = 0;
        std::atomic<float> lastPhoton  = 0;
        parallel_for([&](int index) {
            SPPMPIxel& pixel = pixels[index];
            if (pixel.curPhotonCount > 0) {
                if (iteration > 0) {
                }
                Float factor     = 2.f / 3;
                Float newPhotons = pixel.curPhotonCount * factor + pixel.lastPhotonCount;

                lastPhoton = lastPhoton + pixel.lastPhotonCount;
                curPhoton  = curPhoton + pixel.curPhotonCount;

                Float    newRadius = pixel.radius * sqrt(newPhotons / (pixel.lastPhotonCount + pixel.curPhotonCount));
                Spectrum phi;
                for (int i = 0; i < 3; i++)
                    phi[i] = pixel.phi[i];
                pixel.flux = (pixel.flux + phi * pixel.vp.beta) * (newRadius * newRadius) /
                             (pixel.radius * pixel.radius);
                // pixel.flux = phi;
                pixel.lastPhotonCount = newPhotons;
                pixel.curPhotonCount  = 0;
                pixel.lastRadius      = pixel.radius;
                pixel.radius          = newRadius;

                ++updateCount;
                rate       = rate + pixel.radius / pixel.lastRadius;
                pixel.rate = pixel.radius / pixel.lastRadius;

                for (int i = 0; i < 3; i++)
                    pixel.phi[i] = 0;
            }
            pixel.vp.beta = Spectrum(0);
            if (pixel.vp.event) {
                delete pixel.vp.event->its;
                delete pixel.vp.event;
                pixel.vp.event = nullptr;
            }
        },
                     pixelNum,
                     4096);
        LOGI("Update pixel completed");

        Float minR          = initRadius;
        Float maxR          = 0;
        Float maxPhotons    = 0;
        Float averageRadius = 0;
        for (int i = 0; i < pixelNum; i++) {
            minR       = std::min(pixels[i].radius, minR);
            maxR       = std::max(pixels[i].radius, maxR);
            maxPhotons = std::max(pixels[i].lastPhotonCount, maxPhotons);
            averageRadius += pixels[i].radius / pixelNum;
            if (pixels[i].radius > 0.25) {
            }
        }

        if (iteration + 1 == iterations || (iteration + 1) % writeFrequency == 0) {

            {
                for (int y = 0; y < pixelBounds.y; y++)
                    for (int x = 0; x < pixelBounds.x; x++) {
                        const SPPMPIxel& pixel   = pixels[y * pixelBounds.x + x];
                        int              photonN = photonsPerIteration * (iteration + 1);
                        Spectrum         L(0);
                        L += (pixel.flux) / (photonN * Constant::PI * pixel.radius * pixel.radius);
                        if (iteration == iterations - 1)
                            L += pixel.Ld / (iteration + 1.f);
                        Spectrum c = Spectrum(1 - (pixel.radius - minR) / (maxR - minR));
                        _camera->image->addPixel(x, y, L, true);
                    }
                _camera->image->save(scene.options.outputFileName, 1, false);
            }
        }
        reporter.update(1);
        LOGI("Iteration {0} end", iteration);
    }

    parallel_cleanup();
}

void PhotonMapper::process(const Scene& scene, Sampler& sampler) {
}

PhotonMapper::PhotonMapper(const std::shared_ptr<Camera>& camera, const Json& json) : Integrator(json),
                                                                                      _camera(camera) {
    initRadius          = getOptional(json, "radius", 0.01);
    iterations          = getOptional(json, "interation_num", 16);
    photonsPerIteration = getOptional(json, "photons_per", camera->image->width() * camera->image->height());
    writeFrequency      = getOptional(json, "write_frequency", 128);
}