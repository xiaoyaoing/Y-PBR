#pragma once
#include "../Camera/Camera.hpp"
#include "../Sampler/Sampler.hpp"
#include "../scene.hpp"
#include "../Integrator/PathIntegrator.hpp"
#include "work-queue.hpp"
struct Bucket
    {
        Bucket() : min(0), max(0) { }
        Bucket(const ivec2& min, const ivec2& max) : min(min), max(max) { }

        ivec2 min;
        ivec2 max;
    };


class Render {
public:
    Render(nlohmann::json);
    void sampleImageThread(Sampler * threadSampler);
    void sampleImage();
    void Go();
private:
    std::unique_ptr<Camera> camera ;
    std::unique_ptr<Image>  image;
    std::unique_ptr<Scene>  scene;
    std::unique_ptr<PathIntegrator> integrator;
    std::shared_ptr<Sampler> sampler;
    std::string outputFile;
    size_t thread_num;
    std::unique_ptr<WorkQueue<Bucket>> render_queue;
    const size_t bucketSize =32;
    size_t spp;
};