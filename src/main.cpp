
#include <string>
#include <fstream>
#include "nlohmann/json.hpp"
#include "Camera/Camera.hpp"
#include "Sampler/Sampler.hpp"
#include "scene.hpp"
#include "Integrator/PathIntegrator.hpp"

int main(){
    std::string scene_path="../scenes/standford_box.json";

    std::ifstream scene_file("/Users/yjp/nju/大三下/graphics/Y-PBR/scenes/standford_box.json");

    nlohmann::json j;
    scene_file >> j;
    scene_file.close();


    std::unique_ptr<Camera> camera = std::make_unique<Camera>(j.at("camera").at(0));
    std::unique_ptr<Image>  outputImage = std::make_unique<Image>(j.at("image"));
    std::unique_ptr<Scene>  scene= std::make_unique<Scene>(j);
    std::unique_ptr<PathIntegrator> integrator=std::make_unique<PathIntegrator>(j);

    UniformSampler sampler;

    Ray ray;

    for(size_t y=0;y<outputImage->height;y++)
    {
        for(size_t x=0;x<outputImage->width;x++)
            {
                camera->sampleRay(x,y,outputImage->width,outputImage->height,
                                  ray,sampler.getNext2D());

                outputImage->addPixel(x,y,integrator->integrate(ray,*scene,sampler));
            }
    }


    outputImage->save();

}


