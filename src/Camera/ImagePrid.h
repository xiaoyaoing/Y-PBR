#include "Image.hpp"
#include "Camera.hpp"
#include "IO/FileUtils.hpp"
#pragma once
class ImagePrid{
    std::vector<Image> images;
    int maxPathLength;
public:
    ImagePrid(int pathLength, const Camera & camera): images(pyramidCount(pathLength), Image(camera.image->resoulation(), ToneMap::Pbrt)), maxPathLength(pathLength){
    }
    static inline int pyramidCount(int pathLength)
    {
        return ((pathLength + 1)*(pathLength + 2))/2 - 1;
    }

    static inline int pyramidIndex(int s, int t)
    {
        return pyramidCount(s + t - 2) + t - 1;
    }
    void saveOutPut(const std::string & fileName,Float scale);
    inline  void addPixel(int s,int t,ivec2 pixel,vec3 rgb,bool count = true){
        if(s==2 && t==1 ){
            int k  = 1;
        }
        images[pyramidIndex(s, t)].addPixel(pixel.x, pixel.y, rgb, count);
    }
};