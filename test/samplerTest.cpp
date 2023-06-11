/**
 * @file 
 * @author JunPing Yuan
 * @brief 
 * @version 0.1
 * @date 2023/6/3
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <fstream>
#include "Sampler/SamplerFactory.h"
#include "Camera/Image.hpp"
#include "IO/FileUtils.hpp"

int main() {
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/test/";
    std::string path = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/test/samplerTest.json";
    std::ifstream scene_file(path);
    nlohmann::json j;
    scene_file >> j;
    scene_file.close();

    auto res = getOptional(j, "res", ivec2(512, 512));
    auto spp = getOptional(j, "spp", 1);
    auto sampleNum = getOptional(j, "sample_num", res.x * res.y * spp);

    auto sampler = SamplerFactory::loadSampler(j.at("sampler"), sampleNum, res);

    sampler->startPixel(ivec2(0,0));

    Image image(res);
    for (int i = 0; i < sampleNum; i++) {
        auto sampled = sampler->getNext2D();
        image.addPixel(sampled.x * res.x, sampled.y * res.y, Spectrum(0.2));
        if(!sampler->startNextSample())
            break;
    }
    image.save("res.png", 1);
}