#pragma  once

#include <nlohmann/json.hpp>
#include "../Common/math.hpp"

class Image{



    vec3 getPixel(int x,int y) const ;

    size_t getIndex(size_t x,size_t y) const ;

    std::string  outputFileName;
    std::vector<Float> pixels;
public:
    Image(nlohmann::json json);

    void addPixel(size_t x,size_t y,vec3 rgb)  ;

    void save() const;

    size_t width;
    size_t height;
};