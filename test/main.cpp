
#include "Common/Render.hpp"
#include "iostream"
#include "Texture/BitMapTexture.hpp"

void convert(std::string file){
    BitMapTexture<Spectrum> t(file);
    t.LoadResources();
    ivec2 res(t.width(),t.height());
    Image image(res,ToneMap::LinearOnly);
    for(int i =0;i<t.width();i++)
        for(int j=0;j<t.height();j++){
            auto value = t.getValue(i,j);
            if(luminace(value)>0.2)
                value = Spectrum(1);
            else
                value = Spectrum(0);
            image.addPixel(i,j,value);
        }
    image.save(file,1);
}


int main(int argc, const char *argv[]) {
//    convert("/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/img.png");
//    exit(0);
    spdlog::set_level(spdlog::level::off);
    std::cout << "Rendering Begin";

    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/sssdragon/";
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/test-ball/";
    //FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/bathroom2/";
    //FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/veach-mis/";
    // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/caustic/";
  //    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/classroom/";
   // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/veach-mis/";
    //  FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/hair/";
     FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/curly-hair/";
    // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/water-caustic/";
  //   FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/livingroom/";
    //  FileUtils::WorkingDir =  "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/furball/";
    //  FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/rayTraceOneWeek/";
    //  FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/staircase/";
     // FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/head/";
    //    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/veach-bidir/";
    //FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/cornell-box/";
    FileUtils::WorkingDir = "/Users/yjp/nju/大三下/graphics/offline-render/Y-PBR/example-scenes/teapot/";

    if(argc>0){
        for(int i =1; i<argc;i++)
            Render::renderScene(std::string(argv[i]));
    }
    else {
        Render::renderScene();
    }

}


