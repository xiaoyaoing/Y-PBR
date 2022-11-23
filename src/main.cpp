
#include "IO/FileUtils.hpp"
#include "Common/Render.hpp"
#include "iostream"

int main(){
    spdlog::set_level(spdlog::level::off);

    std::cout<<"lallala";

    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/cornell-box/";
  //  FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/test-ball/";
//    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/veach-mis/";
    //  FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/caustic/";
   // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/classroom/";
    //FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/veach-mis/";
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/hair/";
    std::string s = FileUtils::WorkingDir + "scene.json";
    std::ifstream scene_file(FileUtils::WorkingDir + "scene.json");
//    scene_file.seekg(0, std::ios_base::end);
//    auto size = scene_file.tellg();
//    scene_file.seekg(0);
//    std::string out;
//    out.resize(static_cast<size_t>(size));
//    scene_file.read(&out[0], size);

    nlohmann::json j;
    scene_file >> j;
    scene_file.close();

    Render render(j);
    render.Go();
}


