#include "lodepng/lodepng.h"
#include "Common/math.hpp"


bool saveLdr(const std::string &path, const uint8 *img, int w, int h, int channels);

bool savePng(const std::string &path, const std::vector<unsigned  char> & image , int width, int height, int channels);