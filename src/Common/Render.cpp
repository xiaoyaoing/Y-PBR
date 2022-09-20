#include "Render.hpp"
#include <thread>
#include <spdlog/spdlog.h>
#include "Common/Debug.hpp"
#include "IO/FileUtils.hpp"
Render::Render(nlohmann::json j) {
    camera = std::make_unique<Camera>(j.at("camera"));
    image = std::make_unique<Image>(j.at("camera").at("resolution"));
    scene= std::make_unique<Scene>(j);
    integrator=std::make_unique<PathIntegrator>(j);
    sampler=std::make_shared<UniformSampler>();
    auto renderJson= j["renderer"];
    outputFile = getOptional(renderJson,"output_file",std::string("result"));
    if(outputFile.find(".")!=std::string::npos){
        outputFile = std::string(outputFile.begin(),outputFile.begin()+outputFile.find("."));
    }
    spp = getOptional(renderJson,"spp",32);
}

void Render::sampleImageThread(Sampler * threadSampler){

    Bucket bucket;
    Ray ray;
    Sampler * sampler2 = sampler.get();
    //UniformSampler sampler1 = *(reinterpret_cast<UniformSampler *>(sampler.get()));
  //  std::unique_ptr<Sampler>  threadSampler(sampler->clone());
  //  std::shared_ptr<Sampler> threadSampler = sampler;

    while(render_queue->getWork(bucket)){

        for(size_t y=bucket.min.y;y<bucket.max.y;y++)
        for(size_t x=bucket.min.x;x<bucket.max.x;x++){
            for(uint32 count=0;count<spp;count++)
            {
                camera->sampleRay(x,y,ray,threadSampler->getNext2D());
              //  if(count!=0) image->addPixel(x,y,vec3(0)); else
                Spectrum  pixelResult = integrator->integrate(ray,*scene,*threadSampler);
                image->addPixel(x,y,pixelResult);
            }
            image->dividePixel(x,y,spp);
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
            buckets.emplace_back(Bucket(glm::ivec2(x, y), glm::ivec2(x_end, y_end)));
        }
    }
    
    render_queue=std::make_unique<WorkQueue<Bucket>>(buckets);
    buckets.clear();


    std::function<void(Render*,Sampler * )> f = &Render::sampleImageThread;
    size_t max_threads = std::thread::hardware_concurrency();

    if(DebugConfig::OnlyOneThread) //debug mode  only one  thread
        max_threads = 1;

    spdlog::info("Thread Count {}",max_threads);
    std::vector<std::unique_ptr<std::thread>> threads(max_threads);
    std::vector<std::unique_ptr<Sampler>> samplers(max_threads);
    for(auto & threadSampler : samplers){
        threadSampler = sampler->clone();
    }

    int threadCount=0;
    for (auto& thread : threads)
    {
        thread = std::make_unique<std::thread>(f, this,samplers[threadCount++].get());
    }

    for (auto& thread : threads)
    {
        thread->join();
    }
}
void Render::Go(){
    scene->setUp();
    integrator->Preprocess(*scene,*sampler);
    sampleImage();
    image->postProgress();
    image->savePNG(FileUtils::WorkingDir+outputFile+ std::to_string(spp)+"spp");
    scene->logDebugInfo();
}