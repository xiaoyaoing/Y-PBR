#include "ImageIO.hpp"
#include "FileUtils.hpp"
#include <lodepng/lodepng.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb/stb_image_write.h"

#define TINYEXR_IMPLEMENTATION

#include "tinyexr/tinyexr.h"
#include <iostream>
#include <sstream>
#include <spdlog/spdlog.h>

#if JPEG_AVAILABLE
#include <csetjmp>
#include <jpeglib.h>
#endif

// Image read and write,mainly learned from Tungsten.

static const int GammaCorrection[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
    2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6,
    6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11,
    12, 12, 13, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19,
    19, 20, 21, 21, 22, 22, 23, 23, 24, 25, 25, 26, 27, 27, 28, 29,
    29, 30, 31, 31, 32, 33, 33, 34, 35, 36, 36, 37, 38, 39, 40, 40,
    41, 42, 43, 44, 45, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
    55, 56, 57, 58, 59, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 84, 85, 86, 87, 88, 90,
    91, 92, 93, 95, 96, 97, 99, 100, 101, 103, 104, 105, 107, 108, 109, 111,
    112, 114, 115, 117, 118, 119, 121, 122, 124, 125, 127, 128, 130, 131, 133, 135,
    136, 138, 139, 141, 142, 144, 146, 147, 149, 151, 152, 154, 156, 157, 159, 161,
    162, 164, 166, 168, 169, 171, 173, 175, 176, 178, 180, 182, 184, 186, 187, 189,
    191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221,
    223, 225, 227, 229, 231, 233, 235, 237, 239, 241, 244, 246, 248, 250, 252, 255};

template <typename T>
static T convertToScalar(TexelConversion request, T r, T g, T b, T a, bool haveAlpha)
{
    if (request == TexelConversion::REQUEST_AUTO)
        request = haveAlpha ? TexelConversion::REQUEST_ALPHA : TexelConversion::REQUEST_AVERAGE;

    switch (request)
    {
    case TexelConversion::REQUEST_AVERAGE:
        return (r + g + b) / 3;
    case TexelConversion::REQUEST_RED:
        return r;
    case TexelConversion::REQUEST_GREEN:
        return g;
    case TexelConversion::REQUEST_BLUE:
        return b;
    case TexelConversion::REQUEST_ALPHA:
        return a;
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
static const stbi_io_callbacks istreamCallback = stbi_io_callbacks{
    &stbiReadCallback, &stbiSkipCallback, &stbiEofCallback};

namespace ImageIO
{

    void
    encodeWithState(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height)
    {
        std::vector<unsigned char> png;
        lodepng::State state; // optionally customize this one

        unsigned error = lodepng::encode(png, image, width, height, state);
        if (!error)
            lodepng::save_file(png, filename);

        // if there's an error, display it
        if (error)
            std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }

    void
    encodeOneStep(const char *filename, const std::vector<unsigned char> &image, unsigned width, unsigned height, int channels)
    {
        // Encode the image
        LodePNGColorType types[] = {LCT_GREY, LCT_GREY_ALPHA, LCT_RGB, LCT_RGBA};
        unsigned error = lodepng::encode(filename, image, width, height, types[channels - 1]);

        // if there's an error, display it
        if (error)
            std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }

    static bool savePng(const std::string &path, const std::vector<unsigned char> &image, int width, int height,
                        int channels)
    {
        encodeOneStep(path.c_str(), image, width, height, channels);
        return true;
    }

    bool saveBMP(const std::string &path, const std::vector<unsigned char> &image, int width, int height,
                 int channels)
    {
        // todo
        return false;
    }

    static void nop(void *) {}

    std::unique_ptr<float[]> loadHdr(const std::string &path, TexelConversion request, int &w, int &h)
    {
        std::shared_ptr<std::ifstream> in = std::make_shared<std::ifstream>(FileUtils::getFileFullPath(path), std::ios::binary);
        if (!in)
            return nullptr;

        int channels;
        std::unique_ptr<float[], void (*)(void *)> img(stbi_loadf_from_callbacks(&istreamCallback, in.get(),
                                                                                 &w, &h, &channels, 0),
                                                       stbi_image_free);
        // We only expect Radiance HDR for now, which only has RGB support.
        if (!img || channels != 3)
            return nullptr;
        int targetChannels = (request == TexelConversion::REQUEST_RGB) ? 3 : 1;
        std::unique_ptr<float[]> texels(new float[w * h * targetChannels]);
        if (targetChannels == 3)
        {
            std::memcpy(texels.get(), img.get(), w * h * targetChannels * sizeof(Float));
        }
        return std::move(texels);
    }

    uint8 *loadPng(const std::string &path, int &w, int &h, int &channels)
    {
        int c;
        return stbi_load(path.c_str(), &w, &h, &channels, 4);
        size_t size = FileUtils::getFileSize(path);
        std::shared_ptr<std::ifstream> in = std::make_shared<std::ifstream>(path, std::ios::binary);
        if (!in)
            return nullptr;
        std::unique_ptr<uint8[]> file(new uint8[size_t(size)]);
        FileUtils::streamRead(*in, file.get(), size);
        uint8 *dst;
        uint32 uw, uh;
        lodepng_decode_memory(&dst, &uw, &uh, file.get(), size, LCT_RGBA, 8);
        w = uw;
        h = uh;
        // std::vector<
        //        std::vector < unsigned char > data(w * h * 4);
        //        for(int i=0;i<w * h * 4;i++) data[i] = dst[i];
        //        savePng(FileUtils::getFilePath(FileUtils::WorkingDir+"texture","png",false),data,w,h,4);

        return dst;
    }
#if JPEG_AVAILABLE
    uint8 *loadJpg(const std::string &path, int &w, int &h, int &channels)
    {
        return stbi_load(path.c_str(), &w, &h, &channels, 3);
        struct jpeg_decompress_struct cinfo;
        struct CustomJerr : jpeg_error_mgr
        {
            std::string fileName;
            std::jmp_buf env;
        } jerr;
        jerr.fileName = path;
        if (setjmp(jerr.env))
            return nullptr;

        cinfo.err = jpeg_std_error(&jerr);
        jerr.output_message = [](j_common_ptr cinfo)
        {
            char buffer[JMSG_LENGTH_MAX];
            (*cinfo->err->format_message)(cinfo, buffer);
            spdlog::error("Jpg decoding issue for file '%s': %s\n",
                          static_cast<CustomJerr *>(cinfo->err)->fileName, buffer);
            std::cout.flush();
        };
        jerr.error_exit = [](j_common_ptr cinfo)
        {
            (*cinfo->err->output_message)(cinfo);
            std::longjmp(static_cast<CustomJerr *>(cinfo->err)->env, 1);
        };

        jpeg_create_decompress(&cinfo);

        uint64 fileSize = FileUtils::getFileSize(path);
        std::shared_ptr<std::ifstream> in = std::make_shared<std::ifstream>(path);
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

        uint8 *result(new uint8[cinfo.output_width * cinfo.output_height * 4]);

        std::unique_ptr<unsigned char *[]> scanlines(new unsigned char *[cinfo.output_height]);
        for (unsigned i = 0; i < cinfo.output_height; ++i)
            scanlines[i] = &result[i * cinfo.output_width * 4];
        while (cinfo.output_scanline < cinfo.output_height)
            jpeg_read_scanlines(&cinfo, &scanlines[cinfo.output_scanline], cinfo.output_height - cinfo.output_scanline);

        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);

        for (int y = 0; y < h; ++y)
        {
            for (int x = w - 1; x >= 0; --x)
            {
                result[y * w * 4 + x * 4 + 3] = 0xFF;
                result[y * w * 4 + x * 4 + 2] = result[y * w * 4 + x * 3 + 2];
                result[y * w * 4 + x * 4 + 1] = result[y * w * 4 + x * 3 + 1];
                result[y * w * 4 + x * 4 + 0] = result[y * w * 4 + x * 3 + 0];
            }
        }

        return std::move(result);
    }
#endif

    std::unique_ptr<uint8[]>
    loadLdr(const std::string &path, TexelConversion request, int &w, int &h, bool gammaCorrect)
    {
        auto fullPath = FileUtils::getFileFullPath(path);
        int channels;
        std::unique_ptr<uint8[], void (*)(void *)> img((uint8 *)0, &nop);

        if (path.ends_with("jpg") || path.ends_with("jpeg"))
        {
            img = std::unique_ptr<uint8[], void (*)(void *)>(loadPng(fullPath, w, h, channels), free);
        }

        else if (path.ends_with("png"))
        {
            img = std::unique_ptr<uint8[], void (*)(void *)>(loadPng(fullPath, w, h, channels), free);
        }
        else
        {
            spdlog::error("{} image not supported yet", path);
            return nullptr;
        }
        if (!img)
        {
            spdlog::error("Error when reading {0}", path);
            return nullptr;
        }

        int targetChannels = (request == TexelConversion::REQUEST_RGB) ? 4 : 1;

        std::unique_ptr<uint8[]> texels(new uint8[w * h * targetChannels]);
        if (targetChannels == 4)
        {
            std::memcpy(texels.get(), img.get(), w * h * targetChannels * sizeof(uint8));
            if (gammaCorrect)
                for (int i = 0; i < w * h; ++i)
                    for (int t = 0; t < 3; ++t)
                        texels[i * 4 + t] = GammaCorrection[texels[i * 4 + t]];
        }
        else
        {
            for (int i = 0; i < w * h; ++i)
                texels[i] = convertToScalar(request, int(img[i * 4]), int(img[i * 4 + 1]), int(img[i * 4 + 2]),
                                            int(img[i * 4 + 3]), channels == 4);
        }

        return std::move(texels);
    }
    std::unique_ptr<float[]>
    loadLdrNormalize(const std::string &path, TexelConversion request, int &w, int &h, bool gammaCorrect)
    {
        auto fullPath = FileUtils::getFileFullPath(path);
        int channels;
        std::unique_ptr<uint8[], void (*)(void *)> img((uint8 *)0, &nop);
#if JPEG_AVAILABLE
        if (path.ends_with("jpg") || path.ends_with("jpeg"))
        {
            img = std::unique_ptr<uint8[], void (*)(void *)>(loadJpg(path, w, h, channels), free);
        }
#endif
        if (path.ends_with("png"))
        {
            img = std::unique_ptr<uint8[], void (*)(void *)>(loadPng(fullPath, w, h, channels), free);
        }
        else
        {
            spdlog::error("{} image not supported yet", path);
            return nullptr;
        }
        if (!img)
        {
            spdlog::error("Error when reading {0}", path);
            return nullptr;
        }

        int targetChannels = (request == TexelConversion::REQUEST_RGB) ? 4 : 1;

        std::unique_ptr<uint8[]> texels(new uint8[w * h * targetChannels]);
        std::unique_ptr<float[]> texelsNormalized(new float[w * h * targetChannels]);
        if (targetChannels == 4)
        {
            std::memcpy(texels.get(), img.get(), w * h * targetChannels * sizeof(uint8));
            if (gammaCorrect)
                for (int i = 0; i < w * h; ++i)
                    for (int t = 0; t < 3; ++t)
                        texels[i * 4 + t] = GammaCorrection[texels[i * 4 + t]];
            std::transform(texels.get(), texels.get() + w * h * targetChannels, texelsNormalized.get(),
                           [](uint8_t value)
                           { return static_cast<float>(value) / 255.0f; });
        }
        else
        {
            for (int i = 0; i < w * h; ++i)
                texels[i] = convertToScalar(request, int(img[i * 4]), int(img[i * 4 + 1]), int(img[i * 4 + 2]),
                                            int(img[i * 4 + 3]), channels == 4);
        }

        return std::move(texelsNormalized);
    }

    bool isHdr(const std::string &path)
    {
        auto suffix = FileUtils::getFileSuffix(path);
        return suffix == "hdr" || suffix == "exr";
    }

    bool saveLdr(const std::string &path, const uint8 *img, int w, int h, int channels, bool overwrite)
    {
        auto fullPath = FileUtils::getFileFullPath(path);
        auto suffix = FileUtils::getFileSuffix(path);
        if (suffix.empty())
        {
            return false;
        }
        if (suffix == "png")
        {
            return savePng(FileUtils::getFilePath(fullPath, overwrite), std::vector<uint8_t>(img, img + w * h * channels), w, h, channels);
        }
        return false;
    }

    // See `examples/rgbe2exr/` for more details.
    bool SaveEXR(const float *rgb, int width, int height, const char *outfilename)
    {

        EXRHeader header;
        InitEXRHeader(&header);

        EXRImage image;
        InitEXRImage(&image);

        image.num_channels = 3;

        std::vector<float> images[3];
        images[0].resize(width * height);
        images[1].resize(width * height);
        images[2].resize(width * height);

        // Split RGBRGBRGB... into R, G and B layer
        for (int i = 0; i < width * height; i++)
        {
            images[0][i] = rgb[3 * i + 0];
            images[1][i] = rgb[3 * i + 1];
            images[2][i] = rgb[3 * i + 2];
        }

        float *image_ptr[3];
        image_ptr[0] = &(images[2].at(0)); // B
        image_ptr[1] = &(images[1].at(0)); // G
        image_ptr[2] = &(images[0].at(0)); // R

        image.images = (unsigned char **)image_ptr;
        image.width = width;
        image.height = height;

        header.num_channels = 3;
        header.channels = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);
        // Must be (A)BGR order, since most of EXR viewers expect this channel order.
        strncpy(header.channels[0].name, "B", 255);
        header.channels[0].name[strlen("B")] = '\0';
        strncpy(header.channels[1].name, "G", 255);
        header.channels[1].name[strlen("G")] = '\0';
        strncpy(header.channels[2].name, "R", 255);
        header.channels[2].name[strlen("R")] = '\0';

        header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
        header.requested_pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
        for (int i = 0; i < header.num_channels; i++)
        {
            header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;          // pixel type of input image
            header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel type of output image to be stored in .EXR
        }

        const char *err = NULL; // or nullptr in C++11 or later.
        int ret = SaveEXRImageToFile(&image, &header, outfilename, &err);
        if (ret != TINYEXR_SUCCESS)
        {
            fprintf(stderr, "Save EXR err: %s\n", err);
            FreeEXRErrorMessage(err); // free's buffer for an error message
            return ret;
        }
        printf("Saved exr file. [ %s ] \n", outfilename);

        // free(rgb);

        free(header.channels);
        free(header.pixel_types);
        free(header.requested_pixel_types);
        return true;
    }

    bool saveHdr(const std::string &path, const float *img, int w, int h, int channels, bool overwrite)
    {
        auto fullPath = FileUtils::getFileFullPath(path);
        auto suffix = FileUtils::getFileSuffix(path);
        if (suffix.empty())
        {
            return false;
        }
        if (suffix == "exr")
        {
            fullPath = FileUtils::getFilePath(fullPath, overwrite);
            return SaveEXR(img, w, h, fullPath.c_str());
        }
        return false;
        //        if(suffix == "hdr"){
        //            std::ofstream out(path, std::ios::binary);
        //            if (!out)
        //            {
        //                spdlog::error("Failed to open output file!\n");
        //                return false;
        //            }
        //            stbi_wr
        //            if (!stbi_write_hdr_to_func([](void* context, void* data, int size) -> int {
        //                std::ostream& os = *static_cast<std::ostream*>(context);
        //                os.write(static_cast<const char*>(data), size);
        //                return os.good() ? 1 : 0;
        //            }, &out, w, h, 3, img))
        //            {
        //                std::cerr << "Failed to write HDR data to file!\n";
        //                return;
        //            }
        //        }
    }
}
