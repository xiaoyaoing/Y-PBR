#include "Render.hpp"
#include <thread>
#include <spdlog/spdlog.h>
Render::Render(nlohmann::json j) {
    camera = std::make_unique<Camera>(j.at("camera").at(0));
    image = std::make_unique<Image>(j.at("image"));
    scene= std::make_unique<Scene>(j);
    integrator=std::make_unique<PathIntegrator>(j);
    integrator->Preprocess(*scene,*sampler);
    sampler=std::make_shared <UniformSampler>();
}

void Render::sampleImageThread(){

    Bucket bucket;
    Ray ray;
    std::unique_ptr<Sampler>  threadSampler(sampler->clone());
    while(render_queue->getWork(bucket)){

        for(size_t y=bucket.min.y;y<bucket.max.y;y++)
        for(size_t x=bucket.min.x;x<bucket.max.x;x++){
            for(uint32 count=0;count<camera->sample_count;count++)
            {
                camera->sampleRay(x,y,image->width,image->height,
                                  ray,threadSampler->getNext2D());
                image->addPixel(x,y,integrator->integrate(ray,*scene,*threadSampler));
            }
            image->dividePixel(x,y,camera->sample_count);
        }
    }
}
void Render::sampleImage(){
    std::vector<Bucket> buckets;

    for (size_t x = 0; x < image->width; x += bucketSize)
    {
        size_t x_end = x + bucketSize;
        if (x_end >= image->width) x_end = image->width;
        for (size_t y = 0; y < image->height; y += bucketSize)
        {
            size_t y_end = y + bucketSize;
            if (y_end >= image->height) y_end = image->height;
            buckets.push_back(Bucket(glm::ivec2(x, y), glm::ivec2(x_end, y_end)));
        }
    }
    
    render_queue=std::make_unique<WorkQueue<Bucket>>(buckets);
    buckets.clear();


    std::function<void(Render*)> f = &Render::sampleImageThread;
    size_t max_threads = std::thread::hardware_concurrency();
    spdlog::info("Thread Count {}",max_threads);

    std::vector<std::unique_ptr<std::thread>> threads(max_threads);
    for (auto& thread : threads)
    {
       // std::unique_ptr<Sampler> threadSampler(sampler->clone());
//        std::make_shared<std::thread>(f, this, std::ref(buckets));
        thread = std::make_unique<std::thread>(f, this);
    }

    for (auto& thread : threads)
    {
        thread->join();
    }
}
void Render::Go(){
    sampleImage();
    image->postProgress();
    image->saveTGA();
    image->savePPM();
}