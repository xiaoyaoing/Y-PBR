#define CATCH_CONFIG_MAIN

#include <Camera/Camera.hpp>

#include "catch2/catch.hpp"
#include "iostream"
#include "IO/FileUtils.hpp"
#include "Common/Render.hpp"

TEST_CASE("draw-line-test"){
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/cornell-box/";
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/caustic/";
    std::string s = FileUtils::WorkingDir + "scene.json";
    std::ifstream scene_file(FileUtils::WorkingDir + "scene.json");

    nlohmann::json j;
    scene_file >> j;
    scene_file.close();

  //  Render render(j);
//    for(int i=0;i<100;i++){
//
//    }
}
