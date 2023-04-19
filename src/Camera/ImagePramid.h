#include "Image.hpp"
#include "Camera.hpp"
class ImagePramId{
    std::vector<Image> images;
public:
    ImagePramId(int pathLength,const Camera & camera):images(pyramidCount(pathLength),Image(camera.image->resoulation())){
    }
    static inline int pyramidCount(int pathLength)
    {
        return ((pathLength + 1)*(pathLength + 2))/2 - 1;
    }

    static inline int pyramidIndex(int s, int t)
    {
        return pyramidCount(s + t - 2) + t - 1;
    }
    void saveOutPut(const std::string & fileName){
        for( auto & image:images)
        {   image.postProgress();

            image.savePNG()
        }
    }
    inline  void addPixel(int s,int t,ivec2 pixel,vec3 rgb){

    }
};