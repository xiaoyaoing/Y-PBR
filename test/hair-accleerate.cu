
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
#include "Texture/BitMapTexture.hpp"

#include <tiny-cuda-nn/common_device.h>
#include <tiny-cuda-nn/config.h>
#include <tiny-cuda-nn/encodings/grid.h>
#include <tiny-cuda-nn/encodings/oneblob.h>
#include <concurrent_vector.h>

float img_width;
float img_height;
vec2 img_extent;

template<uint32_t stride>
__global__ void eval_image(uint32_t n_elements, cudaTextureObject_t texture, float *__restrict__ xs_and_ys,
                           float *__restrict__ result) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_elements) return;

    uint32_t output_idx = i * stride;
    uint32_t input_idx = i * 2;

    float4 val = tex2D<float4>(texture, xs_and_ys[input_idx], xs_and_ys[input_idx + 1]);

    result[output_idx + 0] = val.x;
    result[output_idx + 1] = val.y;
    result[output_idx + 2] = val.z;

    for (uint32_t i = 3; i < stride; ++i) {
        result[output_idx + i] = 1;
    }
}

template<typename T>
__global__ void
to_ldr(const uint64_t num_elements, const uint32_t n_channels, const uint32_t stride, const T *__restrict__ in,
       uint8_t *__restrict__ out) {
    const uint64_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= num_elements) return;

    const uint64_t pixel = i / n_channels;
    const uint32_t channel = i - pixel * n_channels;

    out[i] = (uint8_t) (powf(fmaxf(fminf(in[pixel * stride + channel], 1.0f), 0.0f), 1.0f / 2.2f) * 255.0f + 0.5f);
}

template<typename T>
void
save_image(const T *image, int width, int height, int n_channels, int channel_stride, const std::string &filename) {
    tcnn::GPUMemory<uint8_t> image_ldr(width * height * n_channels);
    tcnn::linear_kernel(to_ldr<T>, 0, nullptr, width * height * n_channels, n_channels, channel_stride, image,
                        image_ldr.data());

    std::vector<uint8_t> image_ldr_host(width * height * n_channels);
    CUDA_CHECK_THROW(cudaMemcpy(image_ldr_host.data(), image_ldr.data(), image_ldr.size(), cudaMemcpyDeviceToHost));


    ImageIO::saveLdr(filename.c_str(), image_ldr_host.data(), width, height, n_channels, true);

}

using precision_t = tcnn::network_precision_t;

tcnn::GPUMemory<float> load_image(const std::string &filename, int &width, int &height) {
    // width * height * RGBA
    auto out = ImageIO::loadLdrNormalize(filename.c_str(), TexelConversion::REQUEST_RGB, width, height);
    // float* out = load_stbi(&width, &height, filename.c_str());

    tcnn::GPUMemory<float> result(width * height * 4);
    result.copy_from_host(out.get());
    //free(out); // release memory of image data

    return result;
}


void origin(BitMapTexture<vec3> *BitMaptexture, std::shared_ptr<tcnn::Trainer<float, precision_t, precision_t>> trainer,
            std::shared_ptr<tcnn::NetworkWithInputEncoding<precision_t>> network) {
    int width, height;


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
                                  {"otype", "OneBlob"},
                                  {"n_bins",        32},
                          }},
            {"network",   {
                                  {"otype", "FullyFusedMLP"},
                                  // {"otype", "CutlassMLP"},
                                  {"n_neurons",     64},
                                  {"n_hidden_layers", 4},
                                  {"activation", "ReLU"},
                                  {"output_activation", "None"},
                          }},
    };
    Json encoding_opts = config.value("encoding", Json::object());
    Json loss_opts = config.value("loss", Json::object());
    Json optimizer_opts = config.value("optimizer", Json::object());
    Json network_opts = config.value("network", Json::object());
    std::shared_ptr<tcnn::Loss<precision_t>> loss{tcnn::create_loss<precision_t>(loss_opts)};
    std::shared_ptr<tcnn::Optimizer<precision_t>> optimizer{tcnn::create_optimizer<precision_t>(optimizer_opts)};
    // network = std::make_shared<tcnn::NetworkWithInputEncoding<precision_t>>(2, 3, encoding_opts, network_opts);

    //  trainer = std::make_shared<tcnn::Trainer<float, precision_t, precision_t>>(network, optimizer, loss);

    // Second step: create a cuda texture out of this image. It'll be used to generate training data efficiently on the fly
    tcnn::GPUMemory<float> image = load_image(
            "curly-hair_PT_GROUD_TROUTH.png", width, height);
    cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypePitch2D;
    resDesc.res.pitch2D.devPtr = image.data();
    resDesc.res.pitch2D.desc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
    resDesc.res.pitch2D.width = width;
    resDesc.res.pitch2D.height = height;
    resDesc.res.pitch2D.pitchInBytes = width * 4 * sizeof(float);

    cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.filterMode = cudaFilterModeLinear;
    texDesc.normalizedCoords = true;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.addressMode[1] = cudaAddressModeClamp;
    texDesc.addressMode[2] = cudaAddressModeClamp;

    cudaTextureObject_t texture;
    CUDA_CHECK_THROW(cudaCreateTextureObject(&texture, &resDesc, &texDesc, nullptr));

    uint32_t n_coords = img_width * img_height;
    uint32_t n_coords_padded = tcnn::next_multiple(n_coords, tcnn::BATCH_SIZE_GRANULARITY);

    tcnn::GPUMemory<float> sampled_image(n_coords * 3);
    tcnn::GPUMemory<float> xs_and_ys(n_coords_padded * 2);

    std::vector<float> host_xs_and_ys(n_coords * 2);
    int sampling_height = img_width;
    int sampling_width = img_height;
    for (int y = 0; y < sampling_height; ++y) {
        for (int x = 0; x < sampling_width; ++x) {
            int idx = (y * sampling_width + x) * 2;
            host_xs_and_ys[idx + 0] = (float) (x + 0.5) / (float) sampling_width;
            host_xs_and_ys[idx + 1] = (float) (y + 0.5) / (float) sampling_height;
        }
    }

    xs_and_ys.copy_from_host(host_xs_and_ys.data());



    // Fourth step: train the model by sampling the above image and optimizing an error metric

    // Various constants for the network and optimization
    const uint32_t batch_size = 1 << 18;
    const uint32_t n_training_steps = 10000000;
    const uint32_t n_input_dims = 2; // 2-D image coordinate
    const uint32_t n_output_dims = 3; // RGB color

    cudaStream_t inference_stream;
    CUDA_CHECK_THROW(cudaStreamCreate(&inference_stream));
    cudaStream_t training_stream = inference_stream;

    tcnn::pcg32 rng{1337};

    // Auxiliary matrices for training
    tcnn::GPUMatrix<float> training_target(n_output_dims, batch_size);
    tcnn::GPUMatrix<float> training_batch(n_input_dims, batch_size);

    // Auxiliary matrices for evaluation
    tcnn::GPUMatrix<float> prediction(n_output_dims, n_coords_padded);
    tcnn::GPUMatrix<float> inference_batch(xs_and_ys.data(), n_input_dims, n_coords_padded);


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    float tmp_loss = 0;
    uint32_t tmp_loss_counter = 0;

    std::cout << "Beginning optimization with " << n_training_steps << " training steps." << std::endl;

    uint32_t interval = 10;

    for (uint32_t i = 0; i < 1000; ++i) {
        bool print_loss = i % interval == 0;
        bool visualize_learned_func = i % interval == 0;

        // Compute reference values at random coordinates
        {
            tcnn::generate_random_uniform<float>(training_stream, rng, batch_size * n_input_dims,
                                                 training_batch.data());
            tcnn::linear_kernel(eval_image<n_output_dims>, 0, training_stream, batch_size, texture,
                                training_batch.data(),
                                training_target.data());
        }

        // Training step
        {
            auto ctx = trainer->training_step(training_stream, training_batch, training_target);

            if (i % std::min(interval, (uint32_t) 100) == 0) {
                //tmp_loss += trainer->loss(training_stream, *ctx);
                ++tmp_loss_counter;
            }
        }

        // Debug outputs
        {
            if (print_loss) {
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                std::cout << "Step#" << i << ": " << "loss=" << tmp_loss / (float) tmp_loss_counter << " time="
                          << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]"
                          << std::endl;

                tmp_loss = 0;
                tmp_loss_counter = 0;
            }
            if (visualize_learned_func) {
                network->inference(inference_stream, inference_batch, prediction);
                auto filename = fmt::format("{}.png", i);
                std::cout << "Writing '" << filename << "'... ";
                save_image(prediction.data(), sampling_width, sampling_height, 3, n_output_dims, filename);
                std::cout << "done." << std::endl;
            }
            // Don't count visualizing as part of timing
            // (assumes visualize_learned_pdf is only true when print_loss is true)
            if (print_loss) {
                begin = std::chrono::steady_clock::now();
            }
        }

        if (print_loss && i > 0 && interval < 1000) {
            interval *= 10;
        }
        tcnn::free_all_gpu_memory_arenas();

    }

    // Dump final image if a name was specified
}


Json config = {
        {"loss",      {
                              {"otype",  "RelativeL2"}
                      }},
        {"optimizer", {
                              {"otype",  "Adam"},
                              // {"otype", "Shampoo"},
                              {"learning_rate", 5e-3},
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
                              {"nested", {
                                                 {
                                                         {"n_dims_to_encode", 3},
                                                         {"otype", "HashGrid "}
                                                 }, {
                                                            {"n_dims_to_encode", 3},
                                                            {"otype", "OneBlob"},
                                                            {"n_bins", 32},
                                                    },
                                                 {
                                                         {"n_dims_to_encode", 3},
                                                         {"otype", "OneBlob"},
                                                         {"n_bins", 32},
                                                 }

                                         }}
                      }},
        {"network",   {
                              {"otype",  "FullyFusedMLP"},
                              // {"otype", "CutlassMLP"},
                              {"n_neurons",     64},
                              {"n_hidden_layers", 2},
                              {"activation", "ReLU"},
                              {"output_activation", "None"},
                      }},
};

__global__ void save(vec2 uv, float *__restrict__ result) {
    result[0] = uv.x;
    result[1] = uv.y;

}

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

void cpu_save(vec3 pos, vec3 dir, vec3 tangent, float *result) {
    //  tangent = vec3(0);
    //spdlog::info(dir);
    result[0] = pos[0];
    result[1] = pos[1];
    result[2] = pos[2];
    result[3] = dir[0];
    result[4] = dir[1];
    result[5] = dir[2];
    result[6] = tangent[0];
    result[7] = tangent[1];
    result[8] = tangent[2];

//    result[0] = 0;
//    result[1] = 0;
//    result[2] = 0;
//    result[3] = 0;
//    result[4] = 0;
//    result[5] = 0;
//    result[6] = 0;
//    result[7] = 0;
//    result[8] = 0;

}

void help(cudaStream_t training_stream, uint32_t batch_size, cudaTextureObject_t texture,
          tcnn::GPUMatrix<float> &training_batch,
          tcnn::GPUMatrix<float> &training_target) {
    tcnn::linear_kernel(eval_image<3>, 0, training_stream, batch_size, texture,
                        training_batch.data(),
                        training_target.data());
}

__global__ void save_out(vec3 L, float *__restrict__ result) {
    result[0] = L[0];
    result[1] = L[1];
    result[2] = L[2];
}

__global__ void cuda_copy(uint32_t n_elements, float *__restrict__ src, float *__restrict__ dst) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_elements) return;
    uint32_t idx = 3 * i;
    dst[idx] = src[idx];
    dst[idx + 1] = src[idx + 1];
    dst[idx + 2] = src[idx + 2];
};

const Float uv_scale_factor = 1.f;

class HairIntegrator : public PathIntegrator {
    const uint32_t batch_size = 1 << 18;
    bool by_pos = true;
    bool by_dir = false;
    const uint32_t n_input_dims = 9; //pos,tangent,dir
    const uint32_t n_output_dims = 3;// rgb color
    const int train_num = 512 * 512;
    cudaStream_t training_stream;
    cudaStream_t inference_stream;

    std::shared_ptr<tcnn::NetworkWithInputEncoding<precision_t>> network;

//    GPUMatrix<float> training_target;
//    GPUMatrix<float> training_batch;
    std::atomic<int> train_count;

    std::shared_ptr<tcnn::Trainer<float, precision_t, precision_t>> trainer;
    int beta = 2;
public:
    HairIntegrator(std::shared_ptr<Camera> camera, std::shared_ptr<Sampler> sampler) : PathIntegrator(camera, sampler,
                                                                                                      Json()),
                                                                                       tangent_img(
                                                                                               camera->image->resoulation()),
                                                                                       pos_img(camera->image->resoulation()),
                                                                                       dir_img(camera->image->resoulation()),
                                                                                       hit_hair_img(
                                                                                               camera->image->resoulation()),
                                                                                       hit_hair(
                                                                                               camera->image->resoulation().x *
                                                                                               camera->image->resoulation().y,
                                                                                               false),
                                                                                       training_batch(n_input_dims,
                                                                                                      submit_train_batch),
                                                                                       training_target(n_output_dims,
                                                                                                       submit_train_batch),
                                                                                       predict_batch(n_input_dims,
                                                                                                     submit_prdict_batch),
                                                                                       predict_result(n_output_dims,
                                                                                                      submit_prdict_batch),
                                                                                       predict_pixel_vector(
                                                                                               submit_prdict_batch),
                                                                                       host_predict_result(
                                                                                               submit_prdict_batch *
                                                                                               n_output_dims),
                                                                                       host_predict_value(
                                                                                               submit_prdict_batch *
                                                                                               n_input_dims),
                                                                                       predict_value_vecotr(
                                                                                               submit_prdict_batch *
                                                                                               n_input_dims),
                                                                                       gt("") {
        if (n_input_dims == 2 || n_input_dims == 3) {
            config["encoding"] = {
                    {"otype", "HashGrid"}
            };
        }
        if (n_input_dims == 6) {
            config["encoding"] = {
                    {"nested", {
                            {
                                    {"n_dims_to_encode", 3},
                                    {"otype", "OneBlob"},
                                    {"n_bins", 32400},
                            },
                            {
                                    {"n_dims_to_encode", 3},
                                    {"otype", "OneBlob"},
                                    {"n_bins", 32400},
                            }

                    }}
            };
        }
        if (by_pos && n_input_dims == 3)
            config["encoding"] = {
                    {"otype", "HashGrid"}
            };
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

    Ray getSampleRay(ivec2 res, vec2 u1, vec2 u2) {
        int x = u1.x * res.x;
        int y = u1.y * res.y;
        return _camera->sampleRay(x, y, u2);
    }

    const int submit_train_batch = 256;
    const int submit_prdict_batch = submit_train_batch;
    std::atomic<int> cur_submit_batch = 0;
    std::atomic<int> cur_predict_batch = 0;
    Concurrency::concurrent_vector<vec2> predict_pixel_vector;

    tcnn::GPUMatrix<float> training_batch;
    tcnn::GPUMatrix<float> training_target;

    tcnn::GPUMatrix<float> predict_batch;
    tcnn::GPUMatrix<float> predict_result;
    std::vector<float> host_predict_result;
    std::vector<float> host_predict_value;
    std::vector<bool> hit_hair;
    Concurrency::concurrent_vector<float> predict_value_vecotr;
    BitMapTexture<vec3> gt;
    Image tangent_img, pos_img, dir_img, hit_hair_img;

//    void trainNetWork(const Scene &scene, ivec2 res, Sampler *sampler) {
//        auto ray = getSampleRay(res, sampler);
//        vec3 pos, tangent, dir(0), LPrime(0), L(0);
//
//        ///reutrn pos,dir,tangent,LPrime
//        while (true) {
//            if (integrate(ray, scene, maxBounces, *sampler, pos, dir, tangent, LPrime, L))
//                break;
//            ray = getSampleRay(res, sampler);
//        }
//
//        auto E = L - LPrime;
//        //E = vec3(0);
//        save<<< 1, 1>>>(pos, dir, tangent, training_batch.data() + n_input_dims * cur_submit_batch);
//        save_out<<<1, 1>>>(E, training_target.data() + n_output_dims * cur_submit_batch);
//        cur_submit_batch++;
//        if (cur_submit_batch % submit_train_batch == 0) {
//            cur_submit_batch = 0;
//            trainer->training_step(training_stream, training_batch, training_target);
//        }
//    }

    const int float_size_factor = 4;

    void save_predict() {
        ivec2 res = _camera->image->resoulation();
        std::vector<ivec2> pixel_pos;
        std::vector<float> predict_data;
        for (int i = 0; i < res.x; i++)
            for (int j = 0; j < res.y; j++) {
                if (!isBlack(hit_hair_img.getPixel(i + j * res.x))) {
                    pixel_pos.emplace_back(i, j);
                }
            }
        int hair_pixel_size = (pixel_pos.size() / tcnn::BATCH_SIZE_GRANULARITY) * tcnn::BATCH_SIZE_GRANULARITY;
        predict_data.resize(n_input_dims * hair_pixel_size);
        vec3 min_pos(1e5f);
        vec3 max_pos(-1e5f);
        for (int i = 0; i < hair_pixel_size; i++) {
            auto p = pixel_pos[i];
            int idx = p.x + p.y * res.x;
            if (n_input_dims == 9)
                cpu_save(pos_img.getPixel(idx), dir_img.getPixel(idx), tangent_img.getPixel(idx),
                         predict_data.data() + i * n_input_dims);
            if (n_input_dims == 2) {
                vec2 uv((p.x / img_width) * uv_scale_factor, (1 - p.y / img_height) * uv_scale_factor);
                if (by_pos)
                    uv = pos_img.getPixel(idx);
                (predict_data.data() + i * n_input_dims)[0] = uv.x;
                (predict_data.data() + i * n_input_dims)[1] = uv.y;

            }
            if (n_input_dims == 3) {
                auto vec = tangent_img.getPixel(idx);
                if (by_pos) {
                    vec = pos_img.getPixel(idx);
                    min_pos = min(vec, min_pos);
                    max_pos = max(vec, max_pos);
                }
                if (by_dir) vec = dir_img.getPixel(idx);
                (predict_data.data() + i * n_input_dims)[0] = vec[0];
                (predict_data.data() + i * n_input_dims)[1] = vec[1];
                (predict_data.data() + i * n_input_dims)[2] = vec[2];
            }
            if (n_input_dims == 6) {
                auto tangent = tangent_img.getPixel(idx);
                auto dir = dir_img.getPixel(idx);
                (predict_data.data() + i * n_input_dims)[0] = dir[0];
                (predict_data.data() + i * n_input_dims)[1] = dir[1];
                (predict_data.data() + i * n_input_dims)[2] = dir[2];
                dir = tangent;
                (predict_data.data() + i * n_input_dims)[3] = dir[0];
                (predict_data.data() + i * n_input_dims)[4] = dir[1];
                (predict_data.data() + i * n_input_dims)[5] = dir[2];
            }
        }


        predict_batch = tcnn::GPUMatrix<float>(n_input_dims, hair_pixel_size);
        predict_result = tcnn::GPUMatrix<float>(n_output_dims, hair_pixel_size);
        host_predict_result.resize(n_output_dims * hair_pixel_size);

        CUDA_CHECK_THROW(cudaMemcpy(predict_batch.data(), predict_data.data(), predict_data.size() * float_size_factor,
                                    cudaMemcpyHostToDevice));
        network->inference(inference_stream, predict_batch, predict_result);
        cudaMemcpy(host_predict_result.data(), predict_result.data(), host_predict_result.size() * float_size_factor,
                   cudaMemcpyDeviceToHost);
        // save_image(predict_result.data(), 1024, 1024, 3, n_output_dims, "predict.png");

        //   hit_hair_img.save("hit.png", 1,true);


        for (int i = 0; i < hair_pixel_size; i++) {
            auto pixel = pixel_pos[i];

            auto L = vec3(host_predict_result[3 * i], host_predict_result[3 * i + 1],
                          host_predict_result[3 * i + 2]);
            //    L= vec3(1);
            _camera->image->addPixel(pixel.x, pixel.y, L, false);
        }

    }


    void train_image(const Scene &scene, float &tmp_loss, bool count_loss) {

        // Debug outputs
        ivec2 res = _camera->image->resoulation();
        training_batch = tcnn::GPUMatrix<float>(n_input_dims, train_num);
        training_target = tcnn::GPUMatrix<float>(n_output_dims, train_num);
        std::vector<float> host_traing_batch(train_num * n_input_dims);
        std::vector<float> host_traing_target(train_num * n_output_dims);
        auto sampler = _sampler.get();
        vec3 min_pos(1e5f);
        vec3 max_pos(-1e5f);
        for (int i = 0; i < train_num; i++) {
            int x, y;
            vec2 u1 = sampler->getNext2D(), u2 = sampler->getNext2D();
            auto ray = getSampleRay(res, u1, u2);
            vec3 pos(0), tangent(0), dir(0), LPrime(0), L(0);
            while (true) {
                if (integrate(ray, scene, 1, *sampler, pos, dir, tangent, L)) {
                    //integrate(ray, scene, beta, *sampler, pos, dir, tangent, LPrime);
                    break;
                }
                u1 = sampler->getNext2D(), u2 = sampler->getNext2D();
                ray = getSampleRay(res, u1, u2);
            }
            min_pos = min(min_pos, pos);
            max_pos = max(max_pos, pos);
            vec3 E = (L - LPrime);
            auto uv = (u1 + u2 / img_extent ) * uv_scale_factor;
            if (n_input_dims == 2) {
                //       uv = sampler->getNext2D();
                E = gt.eval(uv);
            }
            E = gt.eval(uv);
             x = uv.x * img_width;
           y = uv.y * img_height;
            if (by_pos)
                uv = pos;

//            _camera->image->addPixel(x,y,E);
        //    pos_img.addPixel(x, y,vec3( uv.x,uv.y,0), true);


            if (n_input_dims == 9)
                cpu_save(pos, dir, tangent, host_traing_batch.data() + n_input_dims * i);
            if (n_input_dims == 2) {
                (host_traing_batch.data() + n_input_dims * i)[0] = uv.x;
                (host_traing_batch.data() + n_input_dims * i)[1] = uv.y;
            }
            if (n_input_dims == 3) {
                auto vec = tangent;
                if (by_pos) vec = pos;
                if (by_dir) vec = dir;
                (host_traing_batch.data() + i * n_input_dims)[0] = vec[0];
                (host_traing_batch.data() + i * n_input_dims)[1] = vec[1];
                (host_traing_batch.data() + i * n_input_dims)[2] = vec[2];
//                (host_traing_batch.data() + n_input_dims * i)[0] = uv.x;
//                (host_traing_batch.data() + n_input_dims * i)[1] = uv.y;
            }
            if (n_input_dims == 6) {
                (host_traing_batch.data() + i * n_input_dims)[0] = dir[0];
                (host_traing_batch.data() + i * n_input_dims)[1] = dir[1];
                (host_traing_batch.data() + i * n_input_dims)[2] = dir[2];
                dir = tangent;
                (host_traing_batch.data() + i * n_input_dims)[3] = dir[0];
                (host_traing_batch.data() + i * n_input_dims)[4] = dir[1];
                (host_traing_batch.data() + i * n_input_dims)[5] = dir[2];
            }
            host_traing_target[n_output_dims * i] = E[0];
            host_traing_target[n_output_dims * i + 1] = E[1];
            host_traing_target[n_output_dims * i + 2] = E[2];
        }
        if ((by_pos && (n_input_dims == 3 || n_input_dims == 2)) || n_input_dims == 9) {
            auto diff_pos = max_pos - min_pos;
            for (int i = 0; i < train_num; i++) {
                host_traing_batch[i * n_input_dims] =
                        (host_traing_batch[i * n_input_dims] - min_pos[0]) / (diff_pos[0]);
                host_traing_batch[i * n_input_dims + 1] =
                        (host_traing_batch[i * n_input_dims + 1] - min_pos[1]) / (diff_pos[1]);
//                host_traing_batch[i * n_input_dims] = uv.x;
//                host_traing_batch[i * n_input_dims + 1] = uv.y;
                if (n_input_dims >=3 )
                    host_traing_batch[i * n_input_dims + 2] =
                            (host_traing_batch[i * n_input_dims + 2] - min_pos[2]) / (diff_pos[2]);
            }
        }
      //  pos_img.linerarNormalize();
   //     pos_img.save("pos.exr", 1, false);
     //   exit(-1);

        CUDA_CHECK_THROW(cudaMemcpy(training_batch.data(), host_traing_batch.data(),
                                    float_size_factor * host_traing_batch.size(), cudaMemcpyHostToDevice));
        auto error = cudaMemcpy(training_target.data(), host_traing_target.data(),
                                float_size_factor * host_traing_target.size(), cudaMemcpyHostToDevice);

        int width = img_width;
        int height = img_height;
//        ImageIO::loadLdrNormalize("curly-hair_PT_GROUD_TROUTH.png", TexelConversion::REQUEST_RGB, width, height);
//        load_image(
//                "curly-hair_PT_GROUD_TROUTH.png", width, height);
//        tcnn::GPUMemory<float> image = load_image(
//                "curly-hair_PT_GROUD_TROUTH.png", width, height);
//        cudaResourceDesc resDesc;
//        memset(&resDesc, 0, sizeof(resDesc));
//        resDesc.resType = cudaResourceTypePitch2D;
//        resDesc.res.pitch2D.devPtr = image.data();
//        resDesc.res.pitch2D.desc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
//        resDesc.res.pitch2D.width = width;
//        resDesc.res.pitch2D.height = height;
//        resDesc.res.pitch2D.pitchInBytes = width * 4 * sizeof(float);
//
//        cudaTextureDesc texDesc;
//        memset(&texDesc, 0, sizeof(texDesc));
//        texDesc.filterMode = cudaFilterModeLinear;
//        texDesc.normalizedCoords = true;
//        texDesc.addressMode[0] = cudaAddressModeClamp;
//        texDesc.addressMode[1] = cudaAddressModeClamp;
//        texDesc.addressMode[2] = cudaAddressModeClamp;

        //cudaTextureObject_t texture;
        //CUDA_CHECK_THROW(cudaCreateTextureObject(&texture, &resDesc, &texDesc, nullptr));
        tcnn::default_rng_t rng{1337};
//          tcnn::generate_random_uniform<float>(training_stream, rng, batch_size * n_input_dims, training_batch.data());
        //  help(training_stream, batch_size, texture, training_batch, training_target);
        std::vector<float> image_host_traing_target(host_traing_target.size(), 1);
        cudaMemcpy(image_host_traing_target.data(), training_target.data(), host_traing_target.size() * 4,
                   cudaMemcpyDeviceToHost);
        auto ctx = trainer->training_step(training_stream, training_batch, training_target);
        if (count_loss)
            tmp_loss += trainer->loss(training_stream, *ctx);

    }

    void initNN() {
        // training_target=    GPUMatrix<float>(n_output_dims, batch_size);
        //  training_batch =    GPUMatrix<float>(n_input_dims, batch_size);

        //   cudaStream_t inference_stream;
        CUDA_CHECK_THROW(cudaStreamCreate(&inference_stream));
        training_stream = inference_stream;


        Json encoding_opts = config.value("encoding", Json::object());
        Json loss_opts = config.value("loss", Json::object());
        Json optimizer_opts = config.value("optimizer", Json::object());
        Json network_opts = config.value("network", Json::object());


        std::shared_ptr<tcnn::Loss<precision_t>> loss{tcnn::create_loss<precision_t>(loss_opts)};
        std::shared_ptr<tcnn::Optimizer<precision_t>> optimizer{tcnn::create_optimizer<precision_t>(optimizer_opts)};
        network = std::make_shared<tcnn::NetworkWithInputEncoding<precision_t>>(n_input_dims, n_output_dims,
                                                                                encoding_opts,
                                                                                network_opts);

        trainer = std::make_shared<tcnn::Trainer<float, precision_t, precision_t>>(network, optimizer, loss);

    }

public:
    void render(const Scene &scene) override {
        // origin(nullptr, trainer, network);
        //  return;
        gt = BitMapTexture<vec3>("curly-hair_PT_GROUD_TROUTH.png");
        gt.LoadResources();
        auto tileSize = scene.options.tileSize;
        ivec2 renderBounds = _camera->image->resoulation();
        int width = _camera->image->width();
        int height = _camera->image->height();
        ivec2 numTiles{(renderBounds.x + tileSize - 1) / tileSize, (renderBounds.y + tileSize - 1) / tileSize};

        int num_threads = std::thread::hardware_concurrency();
        parallel_init(num_threads);

        int spp = scene.options.spp;
        spp = 1000;
        int sppStep = scene.options.sppStep;


        ProgressReporter reporter(numTiles.x * numTiles.y);

        std::vector<std::unique_ptr<Sampler>> samplers(numTiles.x * numTiles.y);
        for (int x = 0; x < numTiles.x; x++)
            for (int y = 0; y < numTiles.y; y++) {
                int seed = y * tileSize + x;
                samplers[seed] = std::move(_sampler->clone(seed));
            }
        /// train one scenond

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        float tmp_loss = 0;
        for (int i = 0; i < spp; i++) {
            int interval = 10;
            bool print_loss = i % interval == 0;
            print_loss = true;
            train_image(scene, tmp_loss, print_loss);

            if (print_loss) {
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                std::cout << "Step#" << i << ": " << "loss=" << tmp_loss << " time="
                          << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]"
                          << std::endl;
                tmp_loss = 0;
            }
            parallel_for([&](const vec2 &tile) {
                //  return ;
                int x0 = tile[0] * tileSize;
                int x1 = std::min(x0 + tileSize, width);
                int y0 = tile[1] * tileSize;
                int y1 = std::min
                        (y0 + tileSize, height);
                auto tileSampler = samplers[tile.y * tileSize + tile.x].get();
                for (int y = y0; y < y1; y++) {
                    for (int x = x0; x < x1; x++) {
                        Ray ray = _camera->sampleRay(x, y, tileSampler->getNext2D());
                        vec3 pos(0), tangent(0), dir(0), LPrime(0), L(0);
                        ///reutrn pos,dir,tangent,LPrime
//                           Spectrum  l = PathIntegrator::integrate(ray,scene,*tileSampler);
//                            _camera->image->addPixel(x, y, l, true);
//                             continue;
                        bool hitHair = true;
                        if (n_input_dims == 9 || n_input_dims == 3 || n_input_dims == 6) {
                            hitHair = integrate(ray, scene, beta, *tileSampler, pos, dir, tangent, L);
                        }
                        if(n_input_dims ==2 )
                            hitHair = integrate(ray, scene,1, *tileSampler, pos, dir, tangent, L);
                        if (hitHair) {
                            //integrate(ray, scene,maxBounces, *tileSampler, pos, dir, tangent, L);
                            tangent_img.addPixel(x, y, tangent, true);
                            if(hasNan(pos)){
                                int k = 1;
                            }
                            pos_img.addPixel(x, y, pos);
                            dir_img.addPixel(x, y, dir);
                            hit_hair_img.addPixel(x, y, Spectrum(1));
                            _camera->image->addPixel(x, y, vec3(0), true);
                        } else _camera->image->addPixel(x, y, L, true);

                    }
                    train_count = 0;
                }

            }, numTiles);
            //  tangent_img.save("curly-hair-tangent.png", 1.f, true);
            auto getFileName = [](std::string name, int i) {
                return name + std::to_string(i) + ".png";
            };
//              tangent_img.save(getFileName("tangent",i) ,1, true);
//                dir_img.save(getFileName("dir",i),1, true);
//              pos_img.save(getFileName("pos",i), 1,true);
//              exit(-1);
            //      pos_img.normalize();
//           dir_img.normalize();
//            tangent_img.normalize();
            pos_img.linerarNormalize();
            save_predict();
            tcnn::free_all_gpu_memory_arenas();
            _camera->image->save(std::to_string(i) + ".png", 1, true);


            pos_img.clear();
            dir_img.clear();
            hit_hair_img.clear();
            pos_img.clear();


        }
        parallel_cleanup();
        _camera->image->save(scene.options.outputFileName, 1.f / spp);

    }

    bool
    integrate(const Ray &ray, const Scene &scene, int maxDepth, Sampler &sampler, vec3 &pos, vec3 &dir, vec3 &tangent,
              vec3 &L) const {
        std::optional<Intersection> its;
        SurfaceEvent surfaceEvent;
        Spectrum thr(1.0);
        bool specularBounce = true;
        Ray _ray(ray);
        int bounces = 0;
        for (bounces = 0;; ++bounces) {

            its = scene.intersect(_ray);

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
//            if (bounces == beta)
//                LPrime = L;
            if (!its.has_value() || bounces >= maxDepth)

                break;


            surfaceEvent = makeLocalScatterEvent(&its.value());
            if (bounces == beta) {
                if (its.has_value()) {
                    tangent = its->tangent.value();
                      tangent = its->Ng;
                    pos = its.value().p;
                    dir = surfaceEvent.wo;
                    dir = ray.d;
                  //  dir = its->Ng;
                }
            }
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
                // if(bounces==1) dir = _ray.d;
                if (russian(bounces, sampler, thr))
                    break;
            }
        }
        if (bounces > 4) {
            int k = 1;
        }
//        if (bounces < beta)
//            LPrime = L;
        //   L = vec3(bounces/10.f);
        return bounces > 0 ;
    }

};


int main(int argc, const char *argv[]) {
    img_width = 300;
    img_height = 1000;
    img_extent = vec2(img_width,img_height);
    FileUtils::WorkingDir = argv[1];
    std::ifstream scene_file(FileUtils::WorkingDir + "scene.json");
    nlohmann::json j;
    scene_file >> j;
    scene_file.close();
    Render render(j);
    render.integrator.reset(new HairIntegrator(render.camera, render.sampler));
    render.Go();
}


//using json = Json;
//int main(int argc, char* argv[]) {
//    try {
//        uint32_t compute_capability = tcnn::cuda_compute_capability();
//        if (compute_capability < tcnn::MIN_GPU_ARCH) {
//            std::cerr
//                    << "Warning: Insufficient compute capability " << compute_capability << " detected. "
//                    << "This program was compiled for >=" << tcnn::MIN_GPU_ARCH << " and may thus behave unexpectedly." << std::endl;
//        }
//
//        if (argc < 2) {
//            std::cout << "USAGE: " << argv[0] << " " << "path-to-image.jpg [path-to-optional-config.json]" << std::endl;
//            std::cout << "Sample EXR files are provided in 'data/images'." << std::endl;
//            return 0;
//        }
//
//        json config = {
//                {"loss", {
//                                 {"otype", "RelativeL2"}
//                         }},
//                {"optimizer", {
//                                 {"otype", "Adam"},
//                                 // {"otype", "Shampoo"},
//                                 {"learning_rate", 1e-2},
//                                 {"beta1", 0.9f},
//                                 {"beta2", 0.99f},
//                                 {"l2_reg", 0.0f},
//                                 // The following parameters are only used when the optimizer is "Shampoo".
//                                 {"beta3", 0.9f},
//                                 {"beta_shampoo", 0.0f},
//                                 {"identity", 0.0001f},
//                                 {"cg_on_momentum", false},
//                                 {"frobenius_normalization", true},
//                         }},
//                {"encoding", {
//                                 {"otype", "OneBlob"},
//                                 {"n_bins", 32},
//                         }},
//                {"network", {
//                                 {"otype", "FullyFusedMLP"},
//                                 // {"otype", "CutlassMLP"},
//                                 {"n_neurons", 64},
//                                 {"n_hidden_layers", 4},
//                                 {"activation", "ReLU"},
//                                 {"output_activation", "None"},
//                         }},
//        };
//
//        if (argc >= 3) {
//            std::cout << "Loading custom json config '" << argv[2] << "'." << std::endl;
//            std::ifstream f{argv[2]};
//            config = json::parse(f, nullptr, true, /*skip_comments=*/true);
//        }
//
//        // First step: load an image that we'd like to learn
//        int width, height;
//        tcnn::GPUMemory<float> image = load_image(argv[1], width, height);
//
//        // Second step: create a cuda texture out of this image. It'll be used to generate training data efficiently on the fly
//        cudaResourceDesc resDesc;
//        memset(&resDesc, 0, sizeof(resDesc));
//        resDesc.resType = cudaResourceTypePitch2D;
//        resDesc.res.pitch2D.devPtr = image.data();
//        resDesc.res.pitch2D.desc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
//        resDesc.res.pitch2D.width = width;
//        resDesc.res.pitch2D.height = height;
//        resDesc.res.pitch2D.pitchInBytes = width * 4 * sizeof(float);
//
//        cudaTextureDesc texDesc;
//        memset(&texDesc, 0, sizeof(texDesc));
//        texDesc.filterMode = cudaFilterModeLinear;
//        texDesc.normalizedCoords = true;
//        texDesc.addressMode[0] = cudaAddressModeClamp;
//        texDesc.addressMode[1] = cudaAddressModeClamp;
//        texDesc.addressMode[2] = cudaAddressModeClamp;
//
//        cudaTextureObject_t texture;
//        CUDA_CHECK_THROW(cudaCreateTextureObject(&texture, &resDesc, &texDesc, nullptr));
//
//        // Third step: sample a reference image to dump to disk. Visual comparison of this reference image and the learned
//        //             function will be eventually possible.
//
//        int sampling_width = width;
//        int sampling_height = height;
//
//        // Uncomment to fix the resolution of the training task independent of input image
//        // int sampling_width = 1024;
//        // int sampling_height = 1024;
//
//        uint32_t n_coords = sampling_width * sampling_height;
//        uint32_t n_coords_padded = tcnn::next_multiple(n_coords, tcnn::batch_size_granularity);
//
//        tcnn::GPUMemory<float> sampled_image(n_coords * 3);
//        tcnn::GPUMemory<float> xs_and_ys(n_coords_padded * 2);
//
//        std::vector<float> host_xs_and_ys(n_coords * 2);
//        for (int y = 0; y < sampling_height; ++y) {
//            for (int x = 0; x < sampling_width; ++x) {
//                int idx = (y * sampling_width + x) * 2;
//                host_xs_and_ys[idx+0] = (float)(x + 0.5) / (float)sampling_width;
//                host_xs_and_ys[idx+1] = (float)(y + 0.5) / (float)sampling_height;
//            }
//        }
//
//        xs_and_ys.copy_from_host(host_xs_and_ys.data());
//
//        tcnn::linear_kernel(eval_image<3>, 0, nullptr, n_coords, texture, xs_and_ys.data(), sampled_image.data());
//
//        save_image(sampled_image.data(), sampling_width, sampling_height, 3, 3, "reference.jpg");
//
//        // Fourth step: train the model by sampling the above image and optimizing an error metric
//
//        // Various constants for the network and optimization
//        const uint32_t batch_size = 1 << 18;
//        const uint32_t n_training_steps = argc >= 4 ? atoi(argv[3]) : 10000000;
//        const uint32_t n_input_dims = 2; // 2-D image coordinate
//        const uint32_t n_output_dims = 3; // RGB color
//
//        cudaStream_t inference_stream;
//        CUDA_CHECK_THROW(cudaStreamCreate(&inference_stream));
//        cudaStream_t training_stream = inference_stream;
//
//        tcnn::default_rng_t rng{1337};
//
//        // Auxiliary matrices for training
//        tcnn::GPUMatrix<float> training_target(n_output_dims, batch_size);
//        tcnn::GPUMatrix<float> training_batch(n_input_dims, batch_size);
//
//        // Auxiliary matrices for evaluation
//        tcnn::GPUMatrix<float> prediction(n_output_dims, n_coords_padded);
//        tcnn::GPUMatrix<float> inference_batch(xs_and_ys.data(), n_input_dims, n_coords_padded);
//
//        json encoding_opts = config.value("encoding", json::object());
//        json loss_opts = config.value("loss", json::object());
//        json optimizer_opts = config.value("optimizer", json::object());
//        json network_opts = config.value("network", json::object());
//
//        std::shared_ptr<tcnn::Loss<precision_t>> loss{tcnn::create_loss<precision_t>(loss_opts)};
//        std::shared_ptr<tcnn::Optimizer<precision_t>> optimizer{tcnn::create_optimizer<precision_t>(optimizer_opts)};
//        std::shared_ptr<tcnn::NetworkWithInputEncoding<precision_t>> network = std::make_shared<tcnn::NetworkWithInputEncoding<precision_t>>(n_input_dims, n_output_dims, encoding_opts, network_opts);
//
//        auto trainer = std::make_shared<tcnn::Trainer<float, precision_t, precision_t>>(network, optimizer, loss);
//
//        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//
//        float tmp_loss = 0;
//        uint32_t tmp_loss_counter = 0;
//
//        std::cout << "Beginning optimization with " << n_training_steps << " training steps." << std::endl;
//
//        uint32_t interval = 10;
//
//        for (uint32_t i = 0; i < n_training_steps; ++i) {
//            bool print_loss = i % interval == 0;
//            bool visualize_learned_func = argc < 5 && i % interval == 0;
//
//            // Compute reference values at random coordinates
//            {
//                tcnn::generate_random_uniform<float>(training_stream, rng, batch_size * n_input_dims, training_batch.data());
//                tcnn::linear_kernel(eval_image<n_output_dims>, 0, training_stream, batch_size, texture, training_batch.data(), training_target.data());
//            }
//
//            // Training step
//            {
//                auto ctx = trainer->training_step(training_stream, training_batch, training_target);
//
//                if (i % std::min(interval, (uint32_t)100) == 0) {
//                    tmp_loss += trainer->loss(training_stream, *ctx);
//                    ++tmp_loss_counter;
//                }
//            }
//
//            // Debug outputs
//            {
//                if (print_loss) {
//                    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//                    std::cout << "Step#" << i << ": " << "loss=" << tmp_loss/(float)tmp_loss_counter << " time=" << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
//
//                    tmp_loss = 0;
//                    tmp_loss_counter = 0;
//                }
//
//                if (visualize_learned_func) {
//                    network->inference(inference_stream, inference_batch, prediction);
//                    auto filename = fmt::format("{}.png", i);
//                    std::cout << "Writing '" << filename << "'... ";
//                    save_image(prediction.data(), sampling_width, sampling_height, 3, n_output_dims, filename);
//                    std::cout << "done." << std::endl;
//                }
//
//                // Don't count visualizing as part of timing
//                // (assumes visualize_learned_pdf is only true when print_loss is true)
//                if (print_loss) {
//                    begin = std::chrono::steady_clock::now();
//                }
//            }
//
//            if (print_loss && i > 0 && interval < 1000) {
//                interval *= 10;
//            }
//        }
//
//        // Dump final image if a name was specified
//        if (argc >= 5) {
//            network->inference(inference_stream, inference_batch, prediction);
//            save_image(prediction.data(), sampling_width, sampling_height, 3, n_output_dims, argv[4]);
//        }
//
//        tcnn::free_all_gpu_memory_arenas();
//
//        // If only the memory arenas pertaining to a single stream are to be freed, use
//        //free_gpu_memory_arena(stream);
//    } catch (std::exception& e) {
//        std::cout << "Uncaught exception: " << e.what() << std::endl;
//    }
//
//    return EXIT_SUCCESS;
//}
//
//
