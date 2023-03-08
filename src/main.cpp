
#include "IO/FileUtils.hpp"
#include "Common/Render.hpp"
#include "iostream"

int main( ) {
    spdlog::set_level(spdlog::level::off);
    std::cout << "Rendering Begin";

    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/cornell-box/";
   // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/sssdragon/";
    //  FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/test-ball/";
      //FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/bathroom2/";
   FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/veach-mis/";
     // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/caustic/";
    // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/classroom/";
    //FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/veach-mis/";
  //  FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/hair/";
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/curly-hair/";
//    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/water-caustic/";
   // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/livingroom/";
 //  FileUtils::WorkingDir =  "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/furball/";


    std::ifstream scene_file(FileUtils::WorkingDir + "scene.json");
    nlohmann::json j;
    scene_file >> j;
    scene_file.close();

    Render render(j);
    render.Go();
}


