#include "ImageIO.hpp"
#include "lodepng/lodepng.h"
#include "iostream"


void encodeWithState(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
    std::vector<unsigned char> png;
    lodepng::State state; //optionally customize this one

    unsigned error = lodepng::encode(png, image, width, height, state);
    if(!error) lodepng::save_file(png, filename);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

void encodeOneStep(const char* filename, const std::vector<unsigned char>& image, unsigned width, unsigned height) {
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

bool savePng(const std::string &path, const std::vector<unsigned  char> & image , int width, int height, int channels){
    encodeOneStep(path.c_str(), image, width, height);
}
