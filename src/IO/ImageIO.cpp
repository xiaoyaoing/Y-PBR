#include "ImageIO.hpp"
#include "lodepng/lodepng.h"
#include "iostream"
#include "fstream"
#include <stbi/stb_image.h>

static int stbiReadCallback(void *user, char *data, int size)
{
    std::istream &in = *static_cast<std::istream *>(user);
    in.read(data, size);
    return int(in.gcount());
}
static void stbiSkipCallback(void *user, int n)
{
    std::istream &in = *static_cast<std::istream *>(user);
    in.seekg(n, std::ios_base::cur);
}
static int stbiEofCallback(void *user)
{
    std::istream &in = *static_cast<std::istream *>(user);
    return in.eof();
}
static const stbi_io_callbacks istreamCallback  = stbi_io_callbacks{
        &stbiReadCallback, &stbiSkipCallback, &stbiEofCallback
};

namespace ImageIO {

    void
    encodeWithState(const char * filename, std::vector < unsigned char > & image, unsigned width, unsigned height) {
        std::vector < unsigned char > png;
        lodepng::State state; //optionally customize this one

        unsigned error = lodepng::encode(png, image, width, height, state);
        if ( ! error ) lodepng::save_file(png, filename);

        //if there's an error, display it
        if ( error ) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }

    void
    encodeOneStep(const char * filename, const std::vector < unsigned char > & image, unsigned width, unsigned height) {
        //Encode the image
        unsigned error = lodepng::encode(filename, image, width, height);

        //if there's an error, display it
        if ( error ) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }

    bool savePng(const std::string & path, const std::vector < unsigned char > & image, int width, int height,
                 int channels) {
        encodeOneStep(path.c_str(), image, width, height);
    }

    std::unique_ptr < float[] > loadHdr(const std::string & path, TexelConversion request, int & w, int & h) {
        std::shared_ptr < std::ifstream > in = std::make_shared < std::ifstream >(path);
        if ( ! in )
            return nullptr;

        int channels;
        std::unique_ptr < float[], void (*)(void *) > img(stbi_loadf_from_callbacks(& istreamCallback, in.get(),
                                                                                    & w, & h, & channels, 0),
                                                          stbi_image_free);

        // We only expect Radiance HDR for now, which only has RGB support.
        if ( ! img || channels != 3 )
            return nullptr;

        int targetChannels = ( request == TexelConversion::REQUEST_RGB ) ? 3 : 1;

        std::unique_ptr < float[] > texels(new float[w * h * targetChannels]);
        if ( targetChannels == 3 ) {
            std::memcpy(texels.get(), img.get(), w * h * targetChannels * sizeof(Float));
        }
//    else {
//        for (int i = 0; i < w*h; ++i)
//            texels[i] = convertToScalar(request, img[i*3], img[i*3 + 1], img[i*3 + 2], 1.0f, false);
//    }

        return std::move(texels);
    }

    std::unique_ptr < uint8[] >
    loadLdr(const std::string & path, TexelConversion request, int & w, int & h, bool gammaCorrect) {
        return std::unique_ptr < uint8[] >();
    }

}
