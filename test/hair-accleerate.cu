
#include "Common/Render.hpp"
#include <iostream>
#include <optional>

#include "Integrator/Integrator.hpp"
#include "Integrator/TraceHelper.h"
#include "SampleRecords/SurfaceScatterEvent.hpp"
#include "Common/ProgressReporter.h"
#include "Common/Parallel.h"
#include "Bsdfs/Reflection.hpp"
#include "Common/Debug.hpp"

#include <tiny-cuda-nn/common_device.h>
#include <tiny-cuda-nn/config.h>
#include <tiny-cuda-nn/encodings/grid.h>
#include <tiny-cuda-nn/encodings/oneblob.h>
#include <concurrent_vector.h>

Json config = {
        {"loss",      {
                              {"otype", "RelativeL2"}
                      }},
        {"optimizer", {
                              {"otype", "Adam"},
                              // {"otype", "Shampoo"},
                              {"learning_rate", 1e-2},
                              {"beta1",           0.9f},
                              {"beta2",      0.99f},
                              {"l2_reg",            0.0f},
                              // The following parameters are only used when the optimizer is "Shampoo".
                              {"beta3", 0.9f},
                              {"beta_shampoo", 0.0f},
                              {"identity", 0.0001f},
                              {"cg_on_momentum", false},
                              {"frobenius_normalization", true},
                      }},
        {"encoding",  {
                              {"otype", "Composite"},

                              {"nested",        {
                                                        {
                                                                {"n_dims_to_encode", 3},
                                                                {"otype", "HashGrid"}
                                                        },
                                                        {
                                                                {"n_dims_to_encode", 3},
                                                                {"otype", "OneBlob"},
                                                                {"n_bins", 32},
                                                        },
                                                        {
                                                                {"n_dims_to_encode", 3},
                                                                {"otype", "OneBlob"},
                                                                {"n_bins", 32},
                                                        }
                                                }},
                      }},
        {"network",   {
                              {"otype", "FullyFusedMLP"},
                              // {"otype", "CutlassMLP"},
                              {"n_neurons",     64},
                              {"n_hidden_layers", 2},
                              {"activation", "ReLU"},
                              {"output_activation", "None"},
                      }},
};

using namespace tcnn;
using precision_t = network_precision_t;

__global__ void save(vec3 pos, vec3 dir, vec3 tangent, float *__restrict__ result) {
    result[0] = pos[0];
    result[1] = pos[1];
    result[2] = pos[2];
    result[3] = dir[0];
    result[4] = dir[1];
    result[5] = dir[2];
    result[6] = tangent[0];
    result[7] = tangent[1];
    result[8] = tangent[2];
}

__global__ void save_out(vec3 L, float *__restrict__ result) {
    result[0] = L[0];
    result[1] = L[1];
    result[2] = L[2];
}

class HairIntegrator : public PathIntegrator {
    const uint32_t batch_size = 1 << 18;
    const uint32_t n_input_dims = 9; //pos,tangent,dir
    const uint32_t n_output_dims = 3;// rgb color
    const int train_num = 256;
    cudaStream_t training_stream;
    cudaStream_t inference_stream;
    std::shared_ptr<NetworkWithInputEncoding<precision_t>> network;

//    GPUMatrix<float> training_target;
//    GPUMatrix<float> training_batch;
    std::atomic<int> train_count;

    std::shared_ptr<Trainer<float, precision_t, precision_t>> trainer;
    int beta = 3;
public:
    HairIntegrator(std::shared_ptr<Camera> camera, std::shared_ptr<Sampler> sampler) : PathIntegrator(camera, sampler,
                                                                                                      Json()),
                                                                                       training_batch(n_input_dims,
                                                                                                      submit_train_batch),
                                                                                       training_target(n_output_dims,
                                                                                                       submit_train_batch),
                                                                                       predict_batch(n_input_dims,submit_prdict_batch),
                                                                                       predict_result(n_output_dims,submit_prdict_batch),
                                                                                       predict_pixel_vector(submit_prdict_batch),
                                                                                       host_predict_result(submit_prdict_batch * n_output_dims)
                                                                                                       {
        spdlog::info("Begin init NN");
        initNN();
        spdlog::info("End init NN");
        maxBounces = 100;
    }

    // void  trainTarget(const Scene & scene,int beta,ivec2 pos,Sampler * sampler ){
    //     while()
    //     auto ray  =  _camera->sampleRay(pos.x,pos.y,sampler->getNext2D());
    //     auto L = integrate(ray,scene,sampler,beta);
    //     auto LPrime = integrate(ray,sceene,sampler,std::numeric_limits<int>::max());
    //     auto E = L - LPrime;
    // }
    Ray getSampleRay(ivec2 res, Sampler *sampler) {
        vec2 samplePos = sampler->getNext2D();
        int x = samplePos.x * res.x;
        int y = samplePos.y * res.y;
        return _camera->sampleRay(x, y, sampler->getNext2D());
    }

    const int submit_train_batch = 128;
    const int submit_prdict_batch = 128 * 16;
    std::atomic<int> cur_submit_batch = 0;
    std::atomic<int> cur_predict_batch = 0;
    Concurrency::concurrent_vector<vec2> predict_pixel_vector;

    GPUMatrix<float> training_batch;
    GPUMatrix<float> training_target;

    GPUMatrix<float> predict_batch;
    GPUMatrix<float> predict_result;
    std::vector<float> host_predict_result;
    void trainNetWork(const Scene &scene, ivec2 res, Sampler *sampler) {
        auto ray = getSampleRay(res, sampler);
        vec3 pos, tangent, dir = ray.d, LPrime(0), L(0);

        ///reutrn pos,dir,tangent,LPrime
        while (true) {
            if (integrate(ray, scene, maxBounces, *sampler, pos, tangent, LPrime, L))
                break;
            ray = getSampleRay(res, sampler);
        }

        auto E = L - LPrime;
        cur_submit_batch++;
        save<<< 1, 1>>>(pos, dir, tangent, training_batch.data() + n_input_dims * cur_submit_batch);
        save_out<<<1, 1>>>(E, training_target.data() + n_output_dims * cur_submit_batch);
        if (cur_submit_batch % submit_train_batch == 0) {
            cur_submit_batch = 0;
            trainer->training_step(training_stream, training_batch, training_target);
        }
    }


    void save_predict() {
        cur_predict_batch = 0;
        network->inference(inference_stream, predict_batch, predict_result);
        CUDA_CHECK_THROW(cudaMemcpy(host_predict_result.data(), predict_result.data(), host_predict_result.size(), cudaMemcpyDeviceToHost));

        for (int i = 0; i < submit_prdict_batch; i++) {
            auto pixel = predict_pixel_vector[i];
            auto L = vec3(host_predict_result[3 * i], host_predict_result[3 * i + 1],
                          host_predict_result[3 * i + 2]);
            _camera->image->addPixel(pixel.x, pixel.y, L, true);
        }
       // predict_pixel_vector.clear();
    }

    void initNN() {
        // training_target=    GPUMatrix<float>(n_output_dims, batch_size);
        //  training_batch =    GPUMatrix<float>(n_input_dims, batch_size);

        //   cudaStream_t inference_stream;
        CUDA_CHECK_THROW(cudaStreamCreate(&inference_stream));
        training_stream = inference_stream;


        json encoding_opts = config.value("encoding", json::object());
        json loss_opts = config.value("loss", json::object());
        json optimizer_opts = config.value("optimizer", json::object());
        json network_opts = config.value("network", json::object());


        std::shared_ptr<Loss<precision_t>> loss{create_loss<precision_t>(loss_opts)};
        std::shared_ptr<Optimizer<precision_t>> optimizer{create_optimizer<precision_t>(optimizer_opts)};
        network = std::make_shared<NetworkWithInputEncoding<precision_t>>(n_input_dims, n_output_dims, encoding_opts,
                                                                          network_opts);
        trainer = std::make_shared<Trainer<float, precision_t, precision_t>>(network, optimizer, loss);

    }

public:
    void render(const Scene &scene) override {
        auto tileSize = scene.options.tileSize;
        ivec2 renderBounds = _camera->image->resoulation();
        int width = _camera->image->width();
        int height = _camera->image->height();
        ivec2 numTiles{(renderBounds.x + tileSize - 1) / tileSize, (renderBounds.y + tileSize - 1) / tileSize};

        int num_threads = std::thread::hardware_concurrency();
        parallel_init(num_threads);

        int spp = scene.options.spp;
        int sppStep = scene.options.sppStep;


        ProgressReporter reporter(numTiles.x * numTiles.y);

        /// train one scenond

        for (int i = 0; i < spp; i++) {
            while (train_count++ < train_num)
                trainNetWork(scene, _camera->image->resoulation(), _sampler.get());
            parallel_for([&](const vec2 &tile) {
                int x0 = tile[0] * tileSize;
                int x1 = std::min(x0 + tileSize, width);
                int y0 = tile[1] * tileSize;
                int y1 = std::min
                        (y0 + tileSize, height);

                int seed = x0 * width + y0;
                std::unique_ptr<Sampler> tileSampler = _sampler->clone(seed);


                for (int y = y0; y < y1; y++) {
                    for (int x = x0; x < x1; x++) {
                        Ray ray = _camera->sampleRay(x, y, tileSampler->getNext2D());
                        vec3 pos, tangent, dir = ray.d, LPrime(0), L(0);
                        ///reutrn pos,dir,tangent,LPrime
                        bool hitHair = integrate(ray, scene, maxBounces, *tileSampler, pos, tangent, LPrime, L);
                        if (hitHair) {
                            save<<<1, 1>>>(pos, dir, tangent, predict_batch.data()+3*cur_predict_batch);
                            if(cur_predict_batch == submit_prdict_batch)
                                int k = 1;
                            predict_pixel_vector[cur_predict_batch] = ivec2(x,y);
                            if(++cur_predict_batch == submit_prdict_batch)
                               save_predict();
                        }
                        _camera->image->addPixel(x, y, L, true);
                    }
                    train_count = 0;
                }

            }, numTiles);
        }
        parallel_cleanup();
        _camera->image->save(scene.options.outputFileName, 1.f / spp);


    }

    bool integrate(const Ray &ray, const Scene &scene, int maxDepth, Sampler &sampler, vec3 &pos, vec3 &tangent,
                   vec3 &LPrime, vec3 &L) const {
        std::optional<Intersection> its;
        SurfaceEvent surfaceEvent;
        Spectrum thr(1.0);
        int bounces = minBounces;
        bool specularBounce = true;
        Ray _ray(ray);
        for (bounces = 0;; ++bounces) {
            if (bounces == beta)
                LPrime = L;
            its = scene.intersect(_ray);
            if (bounces == 0) {
                if (its.has_value()) {
                    tangent = its->tangent.value();
                    pos = its.value().p;
                }
            }
            if (specularBounce && bounces > minBounces) {
                if (its.has_value())
                    L += thr * its->Le(-_ray.d);
                else
                    for (auto light: scene.lights) {
                        if (light->flags == int(LightFlags::Infinite)) {
                            L += thr * light->Le(_ray);
                        }
                    }

            }

            if (!its.has_value() || bounces >= maxDepth)
                break;


            surfaceEvent = makeLocalScatterEvent(&its.value());
            if (its->bsdf->Pure(BSDF_FORWARD)) {
                _ray = surfaceEvent.sctterRay(_ray.d);
            } else {
                if (!its->bsdf->Pure(BSDF_PURE_SPECULR) && bounces < maxDepth - 1) {
                    Spectrum Ld = uniformSampleAllLights
                            (surfaceEvent, scene, sampler, nullptr);  //direct lighting
                    L += thr * Ld;
                }
                surfaceEvent.requestType = BSDF_ALL;
                Spectrum f = its->bsdf->sampleF(surfaceEvent, sampler.getNext2D(), false);
                if (isBlack(f) || surfaceEvent.pdf == 0)
                    break;
                BXDFType flags = surfaceEvent.sampleType;
                specularBounce = (flags & BSDF_SPECULAR) != 0;
                thr *= f / surfaceEvent.pdf;
                _ray = surfaceEvent.sctterRay();
                if (russian(bounces, sampler, thr))
                    break;
            }
        }
        if (bounces > 4) {
            int k = 1;
        }
        return bounces > 0;
    }

};


int main(int argc, const char *argv[]) {
    FileUtils::WorkingDir = argv[1];
    std::ifstream scene_file(FileUtils::WorkingDir + "scene.json");
    nlohmann::json j;
    scene_file >> j;
    scene_file.close();
    Render render(j);
   // render.integrator.reset(new HairIntegrator(render.camera, render.sampler));
    render.Go();
}


