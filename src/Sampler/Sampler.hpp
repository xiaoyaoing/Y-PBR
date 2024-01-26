#pragma once
#include <cstdint>
#include <limits>
#include <memory>

#include "Common/math.hpp"

// sampler interface
class Sampler {

public:
    Sampler() {}

    virtual ~Sampler() = default;

    Sampler(int spp = 4) : samplesPerPixel(spp) {}
    //    uint64_t getSeed() const { return rng.getSeed(); }
    //    void setSeed(uint64_t seed) { rng.setSeed(seed); }

    virtual std::unique_ptr<Sampler> clone(int seed) const = 0;

    virtual Float getNext1D() = 0;

    virtual vec2 getNext2D() = 0;

    virtual void startPixel(const ivec2& point) {
        curPixel  = point;
        curSample = 0;
    }

    virtual bool startNextSample() {
        curSample++;
        return curSample < samplesPerPixel;
    }

    int   samplesPerPixel;
    ivec2 curPixel;
    int   curSample;
};

// uniform distribution sampler

class GlobalSampler : public Sampler {
public:
    // GlobalSampler Public Methods
    //    bool StartNextSample() {
    //
    //    }

    void startPixel(const ivec2& point) override {
        Sampler::startPixel(point);
        dimension           = 0;
        intervalSampleIndex = GetIndexForSample(0);
        // Compute _arrayEndDim_ for dimensions used for array samples
        // arrayStartDim + sampleArray1D.size() + 2 * sampleArray2D.size();
    }

    bool startNextSample() override {
        dimension           = 0;
        intervalSampleIndex = GetIndexForSample(curSample + 1);
        return Sampler::startNextSample();
    }

    // bool SetSampleNumber(int64_t sampleNum);

    Float getNext1D() override {
        if (dimension >= arrayStartDim && dimension < arrayEndDim)
            dimension = arrayEndDim;
        return SampleDimension(intervalSampleIndex, dimension++);
    }

    vec2 getNext2D() override {
        if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim)
            dimension = arrayEndDim;
        vec2 p(SampleDimension(intervalSampleIndex, dimension),
               SampleDimension(intervalSampleIndex, dimension + 1));
        dimension += 2;
        return p;
    }

    GlobalSampler(int samplesPerPixel) : Sampler(samplesPerPixel) {}

    virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;

    virtual Float SampleDimension(int64_t index, int dimension) const = 0;

private:
    // GlobalSampler Private Data
    int              dimension;
    int64_t          intervalSampleIndex;
    static const int arrayStartDim = 5;
    int              arrayEndDim;
};