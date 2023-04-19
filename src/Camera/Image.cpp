#include "Image.hpp"
#include "Common/histogram.hpp"
#include "Common/Json.hpp"
#include "IO/ImageIO.hpp"
#include "IO/FileUtils.hpp"

#include <spdlog/spdlog.h>
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
    return buffers[idx].value() / Float(sampleCounts[idx]);
}

uint32 Image::getIndex(uint32 x, uint32 y) const {
    // std::cout<<3 * x + 3 * width * y<<" ";
    return x + _width * y;
}


void Image::addPixel(uint32 x, uint32 y, vec3 rgb) {
    uint32 idx = getIndex(x, y);
    buffers[idx].add(rgb);
    sampleCounts[idx]++;
}

void Image::dividePixel(uint32 x, uint32 y, uint32 count) {
    //  buffers[getIndex(x,y)].rgb/=Float(count);
}

void Image::postProgress() {

    for (int i = 0; i < buffers.size(); i++) {
        Spectrum rgb = clamp(255.f * ToneMap::toneMap(_tonemapType, buffers[i].value() / Float(sampleCounts[i])),
                             vec3(0), 255.f * vec3(1));
        buffers[i] = rgb;
    }
}

void Image::save(const std::string &fileName) const {
    auto extension = FileUtils::getFileSuffix(fileName);
    if (extension.empty()) {
        spdlog::info("Invalid Path {0}", fileName);
        return;
    } else {
        int byteNum = product() * 3;
        auto isHdr = ImageIO::isHdr(fileName);
        if (isHdr) {
            std::unique_ptr<float[]> hdr(new float[byteNum]);
            for (int i = 0; i < buffers.size(); i++) {
                auto rgb = getPixel(i);
                hdr[3 * i] = rgb.r;
                hdr[3 * i + 1] = rgb.g;
                hdr[3 * i + 2] = rgb.b;
            }
            ImageIO::saveHdr(fileName, hdr.get(), width(), height(), 3);
        } else {
            std::unique_ptr<uint8_t[]> ldr(new uint8_t[byteNum]);
            for (int i = 0; i < buffers.size(); i++) {
                auto rgb = getPixel(i);
                rgb = clamp(255.f * ToneMap::toneMap(_tonemapType, rgb), vec3(0), 255.f * vec3(1));
                ldr[3 * i] = uint8_t(rgb.r);
                ldr[3 * i + 1] = uint8_t(rgb.g);
                ldr[3 * i + 2] = uint8_t(rgb.b);
            }
            ImageIO::saveLdr(fileName, ldr.get(), width(), height(), 3);
        }
    }
}


