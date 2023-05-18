#include "Image.hpp"
#include "Camera.hpp"
#include "IO/FileUtils.hpp"
class ImagePramId{
    std::vector<Image> images;
    int maxPathLength;
public:
    ImagePramId(int pathLength,const Camera & camera):images(pyramidCount(pathLength),Image(camera.image->resoulation())),maxPathLength(pathLength){
    }
    static inline int pyramidCount(int pathLength)
    {
        return ((pathLength + 1)*(pathLength + 2))/2 - 1;
    }

    static inline int pyramidIndex(int s, int t)
    {
        return pyramidCount(s + t - 2) + t - 1;
    }
    void saveOutPut(const std::string & fileName,Float scale){
        auto prefix = FileUtils::getFilePrefix(fileName);
        for (int length = 1; length <= maxPathLength; ++length) {
            for (int t = 1; t <= length + 1; ++t) {
                int s = length + 1 - t;
                int idx = pyramidIndex(s, t);
                char buffer[100];
                sprintf(buffer, "%s-s=%d-t=%d.png", prefix.c_str(),s,t);
                images[idx].save(std::string(buffer),scale,true);
            }
        }
    }
    inline  void addPixel(int s,int t,ivec2 pixel,vec3 rgb){
            if(s==2 && t==1 ){
                int k  = 1;
            }
            images[pyramidIndex(s,t)].addPixel(pixel.x,pixel.y,rgb);
    }
};