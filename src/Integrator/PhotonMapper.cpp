#include "PhotonMapper.hpp"
#include "Common/Parallel.h"
#include "Common/ProgressReporter.h"
#include "Bsdfs/Reflection.hpp"
#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "Sampler/LowDiscrepancy.hpp"
#include "IO/ImageIO.hpp"
#include "IO/FileUtils.hpp"
#include <thread>
//sppm mainly learned from pbrt

struct SPPMPIxel {
    Float radius = 0;
    Float lastRadius = 0;
    Spectrum Ld = Spectrum(0);
    Spectrum flux = Spectrum(0);
    std::atomic < Float > phi[3];
    std::atomic < int > curPhotonCount;
    Float lastPhotonCount;
    bool isDielectric = false;
    Float rate = 0 ;
    Float lastRate = 0;
    ivec2 pos;

    int a,b,c;

    struct VisiblePoint {
        VisiblePoint( ) {}
//        VisiblePoint(const vec3 & p, const vec3 & wo, const BSDF * bsdf, const Spectrum & beta)
//                : p(p), wo(wo), bsdf(bsdf), beta(beta) {}
        SurfaceEvent * event = nullptr ;
        Spectrum throughPut = Spectrum(0);

        ~VisiblePoint( ) {
            //delete event;
        }
    } vp;
};

struct SPPMListNode {
    SPPMPIxel * pixel;
    SPPMListNode * next;

    int length( ) {
        SPPMListNode * node = this;
        int ans = 0;
        while ( node && node->pixel ) {
            ans ++;
            node = node->next;
        }
        return ans;
    }
};


struct AtmoicVec3 {
    std::atomic < vec3 > value;

    AtmoicVec3( ) {
        std::atomic_init(& value, vec3(0));
    }

    void add(const vec3 & v) {
        value = atomic_load(& value) + v;
    }

    void divide(const Float a) {
        value = atomic_load(& value) / a;
    }

    void reset( ) {
        std::atomic_init(& value, vec3(0));

    }
};

static std::atomic < int > visiblePhotonCount;
static std::atomic < int > photonGridNum;
static std::atomic < int > notIntersectNum;
static std::atomic < int > count1;
static std::atomic < int > count2;
static std::atomic < int > count3;
static std::atomic < int > count4;
static std::atomic < int > count5;
static std::atomic < int > count6;
static std::atomic < int > photonInGrid;

vec3 dir;
vec3 photonColor;
vec3 photonRayPos;
std::atomic_int dielectricNum = 0;
int photonDielectricNum = 0;
int photonDielectricBounce0Num = 0;
int dielectribVPNum = 0;
int specularRealVpNum = 0;
AtmoicVec3 averageSPosition;
AtmoicVec3 averageSpecularPixelPosition;

AtmoicVec3 averageSpecularPixelthroughput;
AtmoicVec3 averageSpecularPixelPhi;
std::atomic < Float > averageSpecularPixelPhotons;

AtmoicVec3 photonBounce1AveragePosition;
AtmoicVec3 photonBounce0AveragePosition;
AtmoicVec3 photonBounce0AverageDir;
AtmoicVec3 photonBounce1AverageDir;
AtmoicVec3 averageFlux;

std::atomic < Float > averageSpecularPixelRadius = 0;

int photonSNull = 0;
int photonSD = 0;
int photonDS = 0;
int photonSS = 0;
int photonS = 0;
int photonSRuss = 0;


static vec3 toGridPos(const Bounds3 & bounds, const vec3 & worldPos, const vec3 & gridRes) {
    vec3 clampWorldPos = clamp(worldPos, bounds.pMin, bounds.pMax);
    auto t = ( clampWorldPos - bounds.pMin );
    auto b = bounds.Diagonal();
    auto res = t / b * gridRes;
    return res;
}

//convert 3d grid coordinate to linear index
static inline int hash(const ivec3 & p, int hashSize) {
    return (unsigned int) ( ( p.x * 73856093 ) ^ ( p.y * 19349663 ) ^
                            ( p.z * 83492791 ) ) %
           hashSize;
}

void PhotonMapper::render(const Scene & scene) const {
    ivec2 pixelBounds = _camera->image->resoulation();
    int pixelNum = pixelBounds.x * pixelBounds.y;
    std::unique_ptr < SPPMPIxel[] > pixels(new SPPMPIxel[pixelNum]);
    for ( int i = 0 ; i < pixelNum ; i ++ )
    {pixels[i].radius = initRadius;pixels[i].lastRadius = initRadius;}
    UniformSampler sampler;
    const int tileSize = 16;
    ivec2 nTiles(( pixelBounds.x + tileSize - 1 ) / tileSize,
                 ( pixelBounds.y + tileSize - 1 ) / tileSize);

    int num_threads = std::thread::hardware_concurrency();
    parallel_init(num_threads);
    //generate visible points
    ProgressReporter reporter(iterations);
    for ( int iteration = 0 ; iteration < iterations ; iteration ++ ) {
        parallel_for([&](ivec2 tileIndex) {
            std::unique_ptr < Sampler > tileSampler = sampler.clone();
            tileSampler->setSeed(( tileIndex.y * nTiles.x + tileIndex.x ) * ( iteration + 1 ));
            int x0 = tileIndex.x * tileSize;
            int x1 = std::min(x0 + tileSize - 1, pixelBounds.x - 1);
            int y0 = tileIndex.y * tileSize;
            int y1 = std::min(y0 + tileSize - 1, pixelBounds.y - 1);
            for ( int y = y0 ; y <= y1 ; y ++ )
                for ( int x = x0 ; x <= x1 ; x ++ ) {
                    int pixelIndex = y * pixelBounds.x + x;
                    SPPMPIxel & pixel = pixels[pixelIndex];
                    pixel.pos = ivec2(x, y);
                    if ( iteration > 0 && pixel.radius != 0.25 ) {

                    }
                    Ray ray = _camera->sampleRay(x, y, tileSampler->getNext2D());
                    bool specularBounce = true;
                    Spectrum throughput = Spectrum(1);
                    if ( pixelIndex == 0 )
                        int k = 1;
                    bool isDielectric = false;
                    for ( int bounce = 0 ; bounce < maxBounces ; bounce ++ ) {
                        std::optional < Intersection > its = scene.intersect(ray);
                        if ( ! its ) {
                            for ( auto light: scene.lights );// pixel.Ld += throughput * light->environmentLighting(ray);
                            break;
                        }
                        if ( bounce == 0 && its->bsdf->HasFlag(BSDF_SPECULAR) ) {
                            isDielectric = true;
                        }
                        SurfaceEvent event = makeLocalScatterEvent(& its.value());
                        const BSDF * bsdf = its->bsdf;
                        if ( specularBounce )
                            pixel.Ld += throughput * its->Le(- ray.d);
                        pixel.Ld += throughput * uniformSampleOneLight(event, scene, * tileSampler);
                        bool isDiffuse = bsdf->HasFlag(BSDF_DIFFUSE);
                        bool isGlossy = bsdf->HasFlag(BSDF_GLOSSY);
                        //find visible point
                        if ( isDiffuse || ( isGlossy && bounce == maxBounces - 1 ) ) {
                            pixel.vp.event = new SurfaceEvent(event);
                            pixel.vp.throughPut = throughput;
                            pixel.isDielectric = isDielectric;
                            break;
                        }
                        pixel.vp.event->its->p;
                        if ( bounce < maxBounces - 1 ) {
                            event.requestType = BSDF_ALL;
                            Spectrum f = bsdf->sampleF(event, tileSampler->getNext2D(),false);
                            //if(!specularBounce) f*=AbsCosTheta(event.wi);
                            if ( event.pdf == 0 || isBlack(f) ) break;
                            specularBounce = ( event.sampleType & BSDF_SPECULAR ) != 0;
                            throughput *=   f / event.pdf;
                            if ( luminace(throughput) < 0.25 ) {
                                Float russinanProb = luminace(throughput);
                                if ( tileSampler->getNext1D() < russinanProb ) {
                                    throughput /= luminace(throughput);
                                } else {
                                    break;
                                }
                            }
                            ray = event.sctterRay(event.toWorld(event.wi));
                        }
                    }
                  //  averageSPosition.add(pixel.vp.event->its->p);
                    if ( isDielectric ) {
                        dielectribVPNum ++;
                        if ( pixel.vp.event ) {
                            specularRealVpNum ++;
                          //  averageSPosition.add(pixel.vp.event->its->p);
                        }
                    }
                }

        }, nTiles);

        averageSPosition.divide(specularRealVpNum);
        //grid the visible points
        Bounds3 gridBounds;
        Float maxRadius = 0;
        for ( int i = 0 ; i < pixelNum ; i ++ ) {
             SPPMPIxel & pixel = pixels[i];
            pixel.rate = 0;
            if ( isBlack(pixel.vp.throughPut) )
                continue;
            Bounds3 vpBpunds(pixel.vp.event->its->p);
            vpBpunds = Expand(vpBpunds, pixel.radius);
            gridBounds = Union(gridBounds, vpBpunds);
            maxRadius = std::max(maxRadius, pixel.radius);
        }

        ivec3 gridRes;
        vec3 diag = gridBounds.Diagonal();
        Float maxDiag = maxElement(diag);
        int baseRes = int(maxDiag / maxRadius);
        for ( int i = 0 ; i < 3 ; i ++ )
            gridRes[i] = std::max(int(baseRes * ( diag[i] / maxDiag )), 1);

        std::vector < std::atomic < SPPMListNode *>> grids(pixelNum);
        parallel_for([&](int pixelIndex) {
            SPPMPIxel & pixel = pixels[pixelIndex];
            if ( isBlack(pixel.vp.throughPut) )
                return;
            vec3 pMin = pixel.vp.event->its->p - vec3(pixel.radius);
            vec3 pMax = pixel.vp.event->its->p + vec3(pixel.radius);
            ivec3 gridMin = toGridPos(gridBounds, pMin, gridRes);
            ivec3 gridMax = toGridPos(gridBounds, pMax, gridRes);
            for ( int x = gridMin.x ; x <= gridMax.x ; x ++ )
                for ( int y = gridMin.y ; y <= gridMax.y ; y ++ )
                    for ( int z = gridMin.z ; z <= gridMax.z ; z ++ ) {
                        int gridIndex = hash(ivec3(x, y, z), pixelNum);
                        if ( gridIndex == 437299 ) {
                            toGridPos(gridBounds, pMin, gridRes);
                        }
                        SPPMListNode * node = new SPPMListNode();
                        node->pixel = & pixel;
                        node->next = grids[gridIndex];
                        while ( ! grids[gridIndex].compare_exchange_weak(
                                node->next, node) );
                        photonGridNum ++;
                    }
        }, pixelNum, 4096);

        int count = 0;
        std::vector < int > counts;
        for ( auto & node: grids ) {
            count += atomic_load(& node)->length();
            counts.push_back(atomic_load(& node)->length());
        }


        //trace Photons and accumlate contribution
        std::unique_ptr < Distribution1D > lightPowerDistrib = computeLightPowerDistrib(scene);
        parallel_for([&](int photonIndex) {

            uint64_t haltonIndex =
                    (uint64_t) iteration * (uint64_t) photonsPerIteration +
                    photonIndex + 1;
            int haltonDim = 0;

            bool drawLine = sampler.getNext1D() < 0.00000001;
            drawLine = false;

            Spectrum lineColor = Spectrum(RadicalInverse(haltonDim, haltonIndex),
                                          RadicalInverse(haltonDim + 1, haltonIndex),
                                          RadicalInverse(haltonDim + 2, haltonIndex));
            haltonDim += 3;

            Float lightPdf;
            int lightIndex = lightPowerDistrib->SampleDiscrete(RadicalInverse(haltonDim ++, haltonIndex), & lightPdf);
            auto light = scene.lights[lightIndex];

            vec2 posSample(RadicalInverse(haltonDim, haltonIndex), RadicalInverse(haltonDim + 1, haltonIndex));
            vec2 dirSample(RadicalInverse(haltonDim + 2, haltonIndex), RadicalInverse(haltonDim + 3, haltonIndex));
            haltonDim += 5;
            LightSampleResult lightSample = light->sampleDirect(posSample, dirSample);

            if ( lightSample.lightDirPdf == 0 || lightSample.lightPosPdf == 0 || isBlack(lightSample.radiance) )
                return;
            Ray & photonRay = lightSample.ray;
            dir += photonRay.d / Float(photonsPerIteration);
            photonColor += lightSample.radiance / Float(photonsPerIteration);
            photonRayPos += photonRay.o / Float(photonsPerIteration);

            Spectrum throughput = absDot(photonRay.d, lightSample.lightN) * lightSample.radiance
                                  / ( lightSample.lightDirPdf * lightSample.lightPosPdf );
            bool firstDielectric = false;
            for ( int bounce = 0 ; bounce < 8 ; bounce ++ ) {
                std::optional < Intersection > its = scene.intersect(photonRay);

                //if(photonRay.d.z>0) return ;
                //not hitL
                if ( ! its.has_value() ) {
                    if ( bounce == 0 )
                        notIntersectNum ++;
                    if ( bounce == 1 && firstDielectric )
                        photonSNull ++;
                    break;
                }
                if ( bounce == 0 ) {
                    count1 ++;
                    photonBounce0AveragePosition.add(its->p);
                    photonBounce0AverageDir.add(photonRay.d);
                }

                if ( bounce == 0 && its->bsdf->HasFlag(BSDF_SPECULAR) ) {
                    // photonDielectricBounce0Num ++;
                    firstDielectric = true;
                    photonS ++;
                }
                if ( bounce == 1 ) {
                    if ( its->bsdf->HasFlag(BSDF_SPECULAR) ) {
                        if ( firstDielectric ) photonSS ++;
                        else photonDS ++;
                    } else photonSD ++;
                }
                if ( bounce == 0 ) {
                    drawLine = drawLine & its->bsdf->HasFlag(BSDF_SPECULAR);
                }
                if ( drawLine )
                    _camera->drawLine(photonRay.o, its->p, lineColor);
                //Direct lighting has been considered, so contribution is calculated only when depth is greater than 1
                if ( bounce > 0 ) {
                    if ( bounce == 1 ) {
                        count2 ++;
                        photonBounce1AveragePosition.add(its->p);
                        photonBounce1AverageDir.add(photonRay.d);
                    }
                    if ( firstDielectric ) {

                    }
                    vec3 photonP = its->p;
                    if ( gridBounds.Contains(photonP) ) {
                        photonInGrid ++;
                        vec3 gridPos = toGridPos(gridBounds, photonP, gridRes);
                        int gridIndex = hash(gridPos, pixelNum);
                        for ( SPPMListNode * node = atomic_load(& grids[gridIndex]) ;
                              node != nullptr ; node = node->next ) {
                            if ( iteration > 0 ) {

                            }

                            SPPMPIxel * pixel = node->pixel;
                            SPPMPIxel::VisiblePoint & vp = pixel->vp;
                            SurfaceEvent * event = vp.event;
                            const vec3 & p = event->its->p;
                            if ( length2(p - photonP) > pixel->radius * pixel->radius )
                                continue;
                            event->wi = event->toLocal(- photonRay.d);
                            Spectrum contrib = event->its->bsdf->f(* event, false) * throughput;
                            //f(event->its->bsdf->hasFlag(BSDF_SPECULAR)) contrib /= AbsCosTheta(event->wi);
                            if ( hasNan(contrib) ) {

                            }
                            for ( int i = 0 ; i < 3 ; i ++ )
                                pixel->phi[i] = std::atomic_load(& pixel->phi[i]) + contrib[i];
                            ++ pixel->curPhotonCount;
                            ++ visiblePhotonCount;
                        }
                    } else {

                    }
                }
                const BSDF * photonBsdf = its->bsdf;
                SurfaceEvent event = makeLocalScatterEvent(& its.value());
                event.requestType = BSDF_ALL;
                Spectrum f = photonBsdf->sampleF(event, vec2(RadicalInverse(haltonDim + 0, haltonIndex),
                                                             RadicalInverse(haltonDim + 1, haltonIndex)), true);
                //if(!isSpecualr(event.sampleType)) f *=abs(event.wi.z);
                haltonDim += 2;
                Spectrum throughputNew = throughput * f/ event.pdf;

                Float russProb = std::max((Float) 0, 1 - luminace(throughputNew) / luminace(throughput));
                russProb = 0;
                if ( RadicalInverse(haltonDim ++, haltonIndex) < russProb ) {
                    if ( bounce == 0 ) count3 ++;
                    if ( bounce == 0 && firstDielectric )
                        photonSRuss ++;
                    break;
                }
                throughput = throughputNew / ( 1 - russProb );
                vec3 wi = event.toWorld(event.wi);
                photonRay = event.sctterRay(wi);
                if ( firstDielectric && bounce == 0 ) {

                }
            }
        }, photonsPerIteration, 8192);


        photonBounce1AveragePosition.divide(count2);
        photonBounce0AveragePosition.divide(count1);
        photonBounce1AverageDir.divide(count2);
        photonBounce0AverageDir.divide(count1);

        //update pixel
        std::atomic<float> rate=0;
        std::atomic<int> updateCount = 0;
        std::atomic<float> curPhoton = 0;
        std::atomic<float> lastPhoton = 0;
        parallel_for([&](int index) {
            SPPMPIxel & pixel = pixels[index];
            if ( pixel.curPhotonCount > 0 ) {
                if ( iteration > 0 ) {

                }
                Float factor = 2.f / 3;
                Float newPhotons = pixel.curPhotonCount * factor + pixel.lastPhotonCount;

                lastPhoton= lastPhoton + pixel.lastPhotonCount;
                curPhoton = curPhoton + pixel.curPhotonCount;

                Float newRadius = pixel.radius * sqrt(newPhotons / ( pixel.lastPhotonCount + pixel.curPhotonCount ));
                Spectrum phi;
                for ( int i = 0 ; i < 3 ; i ++ )
                    phi[i] = pixel.phi[i];
                pixel.flux = ( pixel.flux + phi * pixel.vp.throughPut ) * ( newRadius * newRadius ) /
                             ( pixel.radius * pixel.radius );
                int a = pixel.curPhotonCount;
                Float b = pixel.lastPhotonCount;
                pixel.lastPhotonCount = newPhotons;
                pixel.curPhotonCount = 0;
                pixel.lastRadius = pixel.radius;
                pixel.radius = newRadius;

                updateCount++;
                rate= rate + pixel.radius / pixel.lastRadius;
                pixel.rate = pixel.radius / pixel.lastRadius;
                if(abs(pixel.rate - pixel.lastRate)<0.0001){

                }
                pixel.lastRate =pixel.rate;
                pixel.c = pixel.a;
                pixel.a= a;
                pixel.b= b;
                if ( pixel.isDielectric ) {
                    averageFlux.add(pixel.flux);
                    dielectricNum ++;
                    averageSpecularPixelPosition.add(pixel.vp.event->its->p);
                    averageSpecularPixelRadius = atomic_load(& averageSpecularPixelRadius) + pixel.radius;

                    averageSpecularPixelthroughput.add(pixel.vp.throughPut);
                    averageSpecularPixelPhi.add(vec3(atomic_load(& pixel.phi[0]), atomic_load(& pixel.phi[1]),
                                                     atomic_load(& pixel.phi[2])));
                    averageSpecularPixelPhotons = atomic_load(& averageSpecularPixelPhotons) + pixel.lastPhotonCount;
                }
                for ( int i = 0 ; i < 3 ; i ++ )
                    pixel.phi[i] = 0;

            }
            pixel.vp.throughPut = Spectrum(0);
            delete pixel.vp.event;
            pixel.vp.event = nullptr;
        }, pixelNum, 4096);

        Float minR = initRadius;
        Float maxR = 0;
        Float maxPhotons = 0;
        Float averageRadius = 0;
        for ( int i = 0 ; i < pixelNum ; i ++ ) {
            minR = std::min(pixels[i].radius, minR);
            maxR = std::max(pixels[i].radius, maxR);
            maxPhotons = std::max(pixels[i].lastPhotonCount, maxPhotons);
            averageRadius += pixels[i].radius / pixelNum;
            if ( pixels[i].radius > 0.25 ) {

            }
        }


        spdlog::error("{0} {1} {2} {3} {4} {5}",rate/updateCount,visiblePhotonCount,averageRadius,curPhoton,lastPhoton,updateCount);
        visiblePhotonCount = 0;

        averageSpecularPixelPosition.divide(dielectricNum);
        averageSpecularPixelthroughput.divide(dielectricNum);
        averageSpecularPixelPhi.divide(dielectricNum);
        averageSpecularPixelPhotons = atomic_load(& averageSpecularPixelPhotons) / dielectricNum;

        averageSpecularPixelRadius = atomic_load(& averageSpecularPixelRadius) / dielectricNum;
        averageFlux.divide(dielectricNum);
        dielectricNum = 0;
        if ( iteration + 1 == iterations || ( iteration + 1 ) % writeFrequency == 0 ) {

//            spdlog::info(" {0} {1} {2} {3} visiblePhoton {4} {5} {6} {7} {8}",
//                         std::atomic_load(& count3),
//                         std::atomic_load(& notIntersectNum), std::atomic_load(& count2),
//                         std::atomic_load(& photonGridNum), std::atomic_load(& visiblePhotonCount), dielectricNum,
//                         photonDielectricNum, dielectribVPNum, photonDielectricBounce0Num);
            visiblePhotonCount = 0;
            averageFlux.reset();
            visiblePhotonCount = 0;
            dielectricNum = 0;

            // if(iteration == iterations-1)
            {
                for ( int y = 0 ; y < pixelBounds.y ; y ++ )
                    for ( int x = 0 ; x < pixelBounds.x ; x ++ ) {
                        const SPPMPIxel & pixel = pixels[y * pixelBounds.x + x];
                        int photonN = photonsPerIteration * ( iteration + 1 );
                        // photonN = photonsPerIteration;
                        Spectrum L(0);
                        //pixel.Ld / ( iteration + 1.f ); //+
                        L += ( pixel.flux ) / Float(photonN * M_PI * pixel.radius * pixel.radius);
                        if ( iteration == iterations - 1 )
                            L += pixel.Ld / ( iteration + 1.f );
                        Spectrum c = Spectrum(1 - ( pixel.radius - minR ) / ( maxR - minR ));
                        //    c = Spectrum (pixel.flux/Float(photonN));
                        _camera->image->addPixel(x, y,L);
                    }
                _camera->image->postProgress();
                _camera->image->savePNG();
                _camera->image->fill(Spectrum(0));
            }

            bool writeRadius = false;
            if ( writeRadius ) {
                std::vector < unsigned char > image;
                image.resize(4 * pixelBounds.x * pixelBounds.y);

                uint64 offset = 0;
                for ( int i = 0 ; i < pixelNum ; i ++ ) {
                    const SPPMPIxel & pixel = pixels[i];
                    Float c = 1 - ( pixel.radius - minR ) / ( maxR - minR );
                    for ( int bit = 0 ; bit < 3 ; bit ++ )
                        image[offset ++] = c * 255;
                }
                ImageIO::savePng(FileUtils::getFilePath(FileUtils::WorkingDir + "radius", "png", false),
                                 image, pixelBounds.x, pixelBounds.y, 3);
            }
        }
        reporter.update(1);

    }

    parallel_cleanup();

}

void PhotonMapper::process(const Scene & scene, Sampler & sampler) {

}

PhotonMapper::PhotonMapper(const std::shared_ptr < Camera > & camera, const Json & json) : Integrator(json),_camera(camera) {
    initRadius = getOptional(json, "radius", 0.05);
    iterations = getOptional(json, "interation_num", 32);
    photonsPerIteration = getOptional(json, "photons_per",
                                          camera->image->width() * camera->image->height());
    writeFrequency = getOptional(json, "write_frequency", 64);
}

