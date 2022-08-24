#include <string>
#include <fstream>
#include "Common/Render.hpp"
#include <spdlog/spdlog.h>
//int main(){
//    std::string scene_path="../scenes/standford_box.json";
//
//    std::ifstream scene_file("/Users/yjp/nju/大三下/graphics/Y-PBR/scenes/standford_box.json");
//
//    nlohmann::json j;
//    scene_file >> j;
//    scene_file.close();
//
//
//    Render render(j);
//    render.Go();
////    std::unique_ptr<Camera> camera = std::make_unique<Camera>(j.at("camera").at(0));
////    std::unique_ptr<Image>  outputImage = std::make_unique<Image>(j.at("image"));
////    std::unique_ptr<Scene>  scene= std::make_unique<Scene>(j);
////    std::unique_ptr<PathIntegrator> integrator=std::make_unique<PathIntegrator>(j);
////    UniformSampler sampler;
////    integrator->Preprocess(*scene,sampler);
////
////    Ray ray;
////
////    spdlog::info("Sample count per pixel:{0}",camera->sample_count);
////    for(uint32 y=0;y<outputImage->height;y++)
////    {
////        for(uint32 x=0;x<outputImage->width;x++)
////            {
////                for(uint32 count=0;count<camera->sample_count;count++)
////                {
////                camera->sampleRay(x,y,outputImage->width,outputImage->height,
////                                  ray,sampler.getNext2D());
////                outputImage->addPixel(x,y,integrator->integrate(ray,*scene,sampler));
////                }
////                outputImage->dividePixel(x,y,camera->sample_count);
////            }
////    }
////
////    outputImage->postProgress();
////    outputImage->savePPM();
////    outputImage->saveTGA();
//
//
//}
//

