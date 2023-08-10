#include "Image.hpp"
#include "Common/histogram.hpp"
#include "Common/Json.hpp"
#include "IO/ImageIO.hpp"
#include "IO/FileUtils.hpp"


#include <iostream>
#include <fstream>


void Image::saveTXT(const std::string &fileName) const {
    std::ofstream file(fileName + "txt");
    for (unsigned int i = 0; i < _height; ++i) {
        for (unsigned int j = 0; j < _width; ++j) {
            const vec3 rgb = getPixel(j, i);
            file << toColorStr(rgb) << std::endl;
        }
        file.close();
    }
}


vec3 Image::getPixel(int x, int y) const {
    const unsigned int idx = getIndex(x, y);
    return getPixel(idx);
}

vec3 Image::getPixel(int idx) const {
    return sampleCounts[idx] != 0 ? buffers[idx].value()/Float(sampleCounts[idx]) : vec3(0);
}

uint32 Image::getIndex(uint32 x, uint32 y) const {
    // std::cout<<3 * x + 3 * width * y<<" ";
    return x + _width * y;
}


void Image::addPixel(uint32 x, uint32 y, vec3 rgb, bool count) {
    if (x < 0 || x >= width() || y < 0 || y >= height())
        return;
    if(hasNan(rgb)){
        throw("Rdiance NAN");
    }
    if(rgb.r <0 || rgb.g<0 || rgb.z<0){
        //todo
    }
    uint32 idx = getIndex(x, y);
    buffers[idx].add(rgb);
    if (count)
        sampleCounts[idx]++;
}


void Image::save(const std::string &fileName, Float scale, bool overwrite) const {
    auto extension = FileUtils::getFileSuffix(fileName);
    if (extension.empty()) {
        spdlog::info("Invalid Path {0}", fileName);
        return;
    } else {
        auto isHdr = ImageIO::isHdr(fileName);
        if (isHdr) {
            int byteNum = product() * 3;
            std::unique_ptr<float[]> hdr(new float[byteNum]);
            for (int i = 0; i < buffers.size(); i++) {
                auto rgb = getPixel(i) * scale;
               hdr[3 * i] = rgb.r;
                hdr[3 * i + 1] = rgb.g;
                hdr[3 * i + 2] = rgb.b;
            }
            ImageIO::saveHdr( fileName, hdr.get(), width(), height(), 3, overwrite);
        } else {
            int byteNum = product() * 3;
            std::unique_ptr<uint8_t[]> ldr(new uint8_t[byteNum]);
            for (int i = 0; i < buffers.size(); i++) {
                auto rgb = getPixel(i) * scale;
                rgb = clamp(255.f * ToneMap::toneMap(_tonemapType, rgb), vec3(0), 255.f * vec3(1));
                ldr[3 * i] = uint8_t(rgb.r);
                ldr[3 * i + 1] = uint8_t(rgb.g);
                ldr[3 * i + 2] = uint8_t(rgb.b);
            }
            ImageIO::saveLdr(fileName, ldr.get(), width(), height(), 3, overwrite);
        }
    }
}


