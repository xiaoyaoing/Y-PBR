#include <fstream>
#include "Image.hpp"
#include "iostream"

void Image::save() const {
        std::ofstream file(outputFileName);

        file << "P3" << std::endl;
        file << width << " " << height << std::endl;
        file << "255" << std::endl;

        for (unsigned int i = 0; i < height; ++i) {
            for (unsigned int j = 0; j < width; ++j) {
                const vec3 rgb = getPixel(j, i);
                const unsigned int R =
                        std::clamp(static_cast<unsigned int>(255.0f * rgb[0]), 0u, 255u);
                const unsigned int G =
                        std::clamp(static_cast<unsigned int>(255.0f * rgb[1]), 0u, 255u);
                const unsigned int B =
                        std::clamp(static_cast<unsigned int>(255.0f * rgb[2]), 0u, 255u);
                file << R << " " << G << " " << B << std::endl;
            }
        }
        file.close();
}

vec3 Image::getPixel(int x, int y) const {
    const unsigned int idx = getIndex(x,y);
    return vec3 (pixels[idx], pixels[idx + 1], pixels[idx + 2]);
}

size_t Image::getIndex(size_t x, size_t y) const {
   // std::cout<<3 * x + 3 * width * y<<" ";
    return 3 * x + 3 * width * y;
}

void Image::addPixel(size_t x, size_t y, vec3 rgb) {
    auto idx = getIndex(x,y);
    int s=pixels.size();
    pixels[idx]+=rgb.x;
    pixels[idx+1]+=rgb.y;
    pixels[idx+2]+=rgb.z;
}

Image::Image(nlohmann::json j) {
    width=j.at("width");
    height=j.at("height");
    outputFileName=j.at("outputFileName");
    pixels.resize(width*height*3);
}
