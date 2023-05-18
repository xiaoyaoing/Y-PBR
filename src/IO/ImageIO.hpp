#include "lodepng/lodepng.h"
#include "Common/math.hpp"



enum class TexelConversion
{
    REQUEST_RGB,
    REQUEST_AVERAGE,
    REQUEST_RED,
    REQUEST_GREEN,
    REQUEST_BLUE,
    REQUEST_ALPHA,
    REQUEST_AUTO,
};

namespace ImageIO {

    bool saveLdr(const std::string &path, const uint8 *img, int w, int h, int channels, bool overwrite = false);
    bool saveHdr(const std::string &path, const float *img, int w, int h, int channels, bool overwrite = false);
//    static  bool
//    savePng(const std::string & path, const std::vector < unsigned char > & image, int width, int height, int channels);

    std::unique_ptr < float[] > loadHdr(const std::string & path, TexelConversion request, int & w, int & h);

    std::unique_ptr < uint8[] > loadLdr(const std::string & path, TexelConversion request, int & w, int & h,
                                        bool gammaCorrect = true);
    bool isHdr(const std::string & path);
}