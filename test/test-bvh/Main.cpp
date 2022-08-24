#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "Common/Render.hpp"
TEST_CASE("bvh-test"){
     mat4 a(0.f,-1.f,0.f,0.f,
            0.f,0.f,0.f,0.f,
            0.f,0.f,0.f,0.f,
            0.f,0.f,0.f,0.f);
     vec4 b(0.f,4.f,0.f,0.f);


     vec3 c = mult(a,b);
     std::string scene_path="../scenes/standford_box.json";

    std::ifstream scene_file("/Users/yjp/nju/大三下/graphics/Y-PBR/scenes/scene.json");

    nlohmann::json j;
    scene_file >> j;
    scene_file.close();

    Render render(j);
    render.Go();
}


