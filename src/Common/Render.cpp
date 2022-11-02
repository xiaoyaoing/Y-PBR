#include "Render.hpp"

#include "Common/Debug.hpp"
#include "IO/FileUtils.hpp"
#include "Debug.hpp"

#include <Integrator/PhotonMapper.hpp>

#include <thread>
#include <spdlog/spdlog.h>



Render::Render(const Json & json) {
    camera = std::make_unique < Camera >(json.at("camera"));
    scene = std::make_unique < Scene >(json);
    sampler = std::make_shared < UniformSampler >();

    //load integrator
    {
       const Json & integratorJson = json.at("integrator");
       const std::string & type = getOptional(integratorJson,"type",std::string("path_tracer"));
       if(type == "path_tracer")
           integrator = std::make_unique<PathIntegrator>(camera,sampler);
       else if(type == "sppm")
       {
            Float radius = getOptional(integratorJson,"radius",0.05);
            int iterations = getOptional(integratorJson,"interation_num",64);
            int photonsPerIteration = getOptional(integratorJson,"photons_per",camera->image->width() * camera->image->height());
            int writeFrequency = getOptional(integratorJson,"write_frequency",32);
            int maxBounces = getOptional(integratorJson,"max_bounces",8);
            integrator = std::make_unique<PhotonMapper>(camera,iterations,radius,maxBounces,photonsPerIteration,writeFrequency);
       }
    }
}

//void Render::sampleImageThread(Sampler * threadSampler) {
//
//    Bucket bucket;
//    Ray ray;
//    Sampler * sampler2 = sampler.get();
//    //UniformSampler sampler1 = *(reinterpret_cast<UniformSampler *>(sampler.get()));
//    //  std::unique_ptr<Sampler>  threadSampler(sampler->clone());
//    //  std::shared_ptr<Sampler> threadSampler = sampler;
//
//    while ( render_queue->getWork(bucket) ) {
//
//        for ( size_t y = bucket.min.y ; y < bucket.max.y ; y ++ )
//            for ( size_t x = bucket.min.x ; x < bucket.max.x ; x ++ ) {
//                if ( DebugConfig::sampleInRange ) {
//                    if ( x < DebugConfig::sampleMinx || x > DebugConfig::sampleMaxx ||
//                         y < DebugConfig::sampleMiny || y > DebugConfig::sampleMaxy ) {
//                        image->addPixel(x, y, vec3(1));
//                        continue;
//                    }
//                }
//                for ( uint32 count = 0 ; count < spp ; count ++ ) {
//                    camera->sampleRay(x, y, ray, threadSampler->getNext2D());
//                    //  if(count!=0) image->addPixel(x,y,vec3(0)); else
//                    Spectrum pixelResult = integrator->integrate(ray, * scene, * threadSampler);
//                    image->addPixel(x, y, pixelResult);
//                }
//                image->dividePixel(x, y, spp);
//            }
//    }
//}
//
//void Render::sampleImage( ) {
//    std::vector < Bucket > buckets;
//
//    for ( size_t x = 0 ; x < image->_width ; x += bucketSize ) {
//        size_t x_end = x + bucketSize;
//        if ( x_end >= image->_width ) x_end = image->_width;
//        for ( size_t y = 0 ; y < image->_height ; y += bucketSize ) {
//            size_t y_end = y + bucketSize;
//            if ( y_end >= image->_height ) y_end = image->_height;
//            buckets.emplace_back(Bucket(glm::ivec2(x, y), glm::ivec2(x_end, y_end)));
//        }
//    }
//
//    render_queue = std::make_unique < WorkQueue < Bucket>>(buckets);
//    buckets.clear();
//
//
//    std::function < void(Render *, Sampler *) > f = & Render::sampleImageThread;
//    size_t max_threads = std::thread::hardware_concurrency();
//
//    if ( DebugConfig::OnlyOneThread ) //debug mode  only one  thread
//        max_threads = 1;
//
//    spdlog::info("Thread Count {}", max_threads);
//    std::vector < std::unique_ptr < std::thread>> threads(max_threads);
//    std::vector < std::unique_ptr < Sampler>> samplers(max_threads);
//    for ( auto & threadSampler: samplers ) {
//        threadSampler = sampler->clone();
//    }
//
//    int threadCount = 0;
//    for ( auto & thread: threads ) {
//        thread = std::make_unique < std::thread >(f, this, samplers[threadCount ++].get());
//    }
//
//    for ( auto & thread: threads ) {
//        thread->join();
//    }
//}

void Render::Go( ) {
    scene->build();
    integrator->process(* scene, * sampler);
    integrator->render(*scene);
}