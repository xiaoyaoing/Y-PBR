#include "ImageIO.hpp"
#include "FileUtils.hpp"
#include <lodepng/lodepng.h>
#include <stbi/stb_image.h>

#include <iostream>
#include <sstream>
#include <spdlog/spdlog.h>

#if JPEG_AVAILABLE
#include <csetjmp>
#include <jpeglib.h>
#endif

//Image read and write,mainly learned from Tungsten.

static const int GammaCorrection[] = {
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   2,   2,   2,
        2,   2,   3,   3,   3,   3,   3,   4,   4,   4,   4,   5,   5,   5,   5,   6,
        6,   6,   7,   7,   7,   8,   8,   8,   9,   9,   9,  10,  10,  10,  11,  11,
        12,  12,  13,  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18,  19,
        19,  20,  21,  21,  22,  22,  23,  23,  24,  25,  25,  26,  27,  27,  28,  29,
        29,  30,  31,  31,  32,  33,  33,  34,  35,  36,  36,  37,  38,  39,  40,  40,
        41,  42,  43,  44,  45,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,
        55,  56,  57,  58,  59,  60,  61,  62,  63,  65,  66,  67,  68,  69,  70,  71,
        72,  73,  74,  75,  77,  78,  79,  80,  81,  82,  84,  85,  86,  87,  88,  90,
        91,  92,  93,  95,  96,  97,  99, 100, 101, 103, 104, 105, 107, 108, 109, 111,
        112, 114, 115, 117, 118, 119, 121, 122, 124, 125, 127, 128, 130, 131, 133, 135,
        136, 138, 139, 141, 142, 144, 146, 147, 149, 151, 152, 154, 156, 157, 159, 161,
        162, 164, 166, 168, 169, 171, 173, 175, 176, 178, 180, 182, 184, 186, 187, 189,
        191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221,
        223, 225, 227, 229, 231, 233, 235, 237, 239, 241, 244, 246, 248, 250, 252, 255
};


template<typename T>
static T convertToScalar(TexelConversion request, T r, T g, T b, T a, bool haveAlpha)
{
    if (request == TexelConversion::REQUEST_AUTO)
        request = haveAlpha ? TexelConversion::REQUEST_ALPHA : TexelConversion::REQUEST_AVERAGE;

    switch(request) {
        case TexelConversion::REQUEST_AVERAGE: return (r + g + b)/3;
        case TexelConversion::REQUEST_RED:     return r;
        case TexelConversion::REQUEST_GREEN:   return g;
        case TexelConversion::REQUEST_BLUE:    return b;
        case TexelConversion::REQUEST_ALPHA:   return a;
        default:
            return T(0);
    }
}

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

    bool saveBMP(const std::string & path, const std::vector < unsigned char > & image, int width, int height,
                 int channels) {
        //todo
    }

    static void nop(void *){}

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
        return std::move(texels);
    }

    uint8 * loadPng(const std::string & path,int & w,int & h,int & channels){
        size_t size = FileUtils::getFileSize(path);
        std::shared_ptr < std::ifstream > in = std::make_shared < std::ifstream >(path);
        if(!in) return nullptr;
        std::unique_ptr<uint8[]> file(new uint8[size_t(size)]);
        FileUtils::streamRead(*in,file.get(),size);
        uint8 * dst;
        uint32 uw,uh;
        lodepng_decode_memory(&dst,&uw,&uh,file.get(),size, LCT_RGBA, 8);
        return dst;
    }
#if JPEG_AVAILABLE
uint8 * loadJpg(const std::string & path,int & w,int & h,int & channels){

        struct jpeg_decompress_struct cinfo;
        struct CustomJerr : jpeg_error_mgr {
            std::string fileName;
            std::jmp_buf env;
        } jerr;
        jerr.fileName = path;
        if (setjmp(jerr.env))
            return nullptr;

        cinfo.err = jpeg_std_error(&jerr);
        jerr.output_message = [](j_common_ptr cinfo) {
            char buffer[JMSG_LENGTH_MAX];
            (*cinfo->err->format_message) (cinfo, buffer);
            spdlog::error("Jpg decoding issue for file '%s': %s\n",
                                     static_cast<CustomJerr *>(cinfo->err)->fileName, buffer);
            std::cout.flush();
        };
        jerr.error_exit = [](j_common_ptr cinfo) {
            (*cinfo->err->output_message)(cinfo);
            std::longjmp(static_cast<CustomJerr *>(cinfo->err)->env, 1);
        };

        jpeg_create_decompress(&cinfo);

        uint64 fileSize = FileUtils::getFileSize(path);
        std::shared_ptr < std::ifstream > in = std::make_shared < std::ifstream >(path);
        if (!in)
            return nullptr;

        // Could use more sophisticated JPEG input handler here, but they are a pain to implement
        // For now, just read the file straight to memory
        std::unique_ptr<uint8[]> fileBuf(new uint8[fileSize]);
        FileUtils::streamRead(*in, fileBuf.get(), fileSize);
        jpeg_mem_src(&cinfo, fileBuf.get(), fileSize);

        jpeg_read_header(&cinfo, TRUE);
        cinfo.out_color_space = JCS_RGB;
        jpeg_start_decompress(&cinfo);

        w = cinfo.output_width;
        h = cinfo.output_height;
        channels = 4;

        cinfo.out_color_space = JCS_RGB;

        uint8 * result(new uint8[cinfo.output_width*cinfo.output_height*4]);

        std::unique_ptr<unsigned char *[]> scanlines(new unsigned char *[cinfo.output_height]);
        for (unsigned i = 0; i < cinfo.output_height; ++i)
            scanlines[i] = &result[i*cinfo.output_width*4];
        while (cinfo.output_scanline < cinfo.output_height)
            jpeg_read_scanlines(&cinfo, &scanlines[cinfo.output_scanline], cinfo.output_height - cinfo.output_scanline);

        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);

        for (int y = 0; y < h; ++y) {
            for (int x = w - 1; x >= 0; --x) {
                result[y*w*4 + x*4 + 3] = 0xFF;
                result[y*w*4 + x*4 + 2] = result[y*w*4 + x*3 + 2];
                result[y*w*4 + x*4 + 1] = result[y*w*4 + x*3 + 1];
                result[y*w*4 + x*4 + 0] = result[y*w*4 + x*3 + 0];
            }
        }

        return std::move(result);
    }
#endif

    std::unique_ptr < uint8[] >
    loadLdr(const std::string & path, TexelConversion request, int & w, int & h, bool gammaCorrect) {
        int channels;
        std::unique_ptr<uint8[], void(*)(void *)> img((uint8 *)0,&nop);
#if JPEG_AVAILABLE
        if(path.ends_with("jpg") || path.ends_with("jpeg")){
            img = std::unique_ptr<uint8[], void(*)(void *)>(loadJpg(path,w,h,channels),free);
        }
#endif
        if(path.ends_with("png")){
            img = std::unique_ptr<uint8[], void(*)(void *)>(loadPng(path,w,h,channels),free);
        }
        else {
            spdlog::error("{} image not supported yet",path);
            return nullptr;
        }
        if(!img){
            spdlog::error("Error when reading {0}",path);
            return nullptr;
        }

        int targetChannels = (request == TexelConversion::REQUEST_RGB) ? 4 : 1;

        std::unique_ptr<uint8[]> texels(new uint8[w*h*targetChannels]);
        if (targetChannels == 4) {
            std::memcpy(texels.get(), img.get(), w*h*targetChannels*sizeof(uint8));
            if (gammaCorrect)
                for (int i = 0; i < w*h; ++i)
                    for (int t = 0; t < 3; ++t)
                        texels[i*4 + t] = GammaCorrection[texels[i*4 + t]];
        } else {
            for (int i = 0; i < w*h; ++i)
                texels[i] = convertToScalar(request, int(img[i*4]), int(img[i*4 + 1]), int(img[i*4 + 2]),
                                            int(img[i*4 + 3]), channels == 4);
        }

        return std::move(texels);
    }

}