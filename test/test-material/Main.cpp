#define CATCH_CONFIG_MAIN

#include "IO/FileUtils.hpp"
#include "catch2/catch.hpp"
#include "Common/Render.hpp"
TEST_CASE("bvh-test"){

    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/Y-PBR/example-scenes/cornell-box/";
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/Y-PBR/example-scenes/test-ball/";

    std::ifstream scene_file(FileUtils::WorkingDir + "scene.json");

     nlohmann::json j;
     scene_file >> j;
     scene_file.close();

     Render render(j);
     render.Go();
}


