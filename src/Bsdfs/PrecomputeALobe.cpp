//
// Created by pc on 2023/7/21.
//

#include "PrecomputeALobe.h"

PrecomputedAzimuthalLobe::PrecomputedAzimuthalLobe(std::unique_ptr<vec3[]> table)
        : _table(std::move(table))
{
    const int Size = AzimuthalResolution;

    std::vector<float> weights(Size*Size);
    for (int i = 0; i < Size*Size; ++i)
        weights[i] = max(_table[i]);

    // Dilate weights slightly to stay conservative
    for (int y = 0; y < Size; ++y) {
        for (int x = 0; x < Size - 1; ++x)
            weights[x + y*Size] = std::max(weights[x + y*Size], weights[x + 1 + y*Size]);
        for (int x = Size - 1; x > 0; --x)
            weights[x + y*Size] = std::max(weights[x + y*Size], weights[x - 1 + y*Size]);
    }
    for (int x = 0; x < Size; ++x) {
        for (int y = 0; y < Size - 1; ++y)
            weights[x + y*Size] = std::max(weights[x + y*Size], weights[x + (y + 1)*Size]);
        for (int y = Size - 1; y > 0; --y)
            weights[x + y*Size] = std::max(weights[x + y*Size], weights[x + (y - 1)*Size]);
    }

    _sampler.reset(new InterpolatedDistribution1D(std::move(weights), Size, Size));
}

void PrecomputedAzimuthalLobe::sample(float cosThetaD, float xi, float &phi, float &pdf) const {


        float v = (AzimuthalResolution - 1) * cosThetaD;

        int x;
        _sampler->warp(v, xi, x);

        phi = Constant::TWO_PI * (x + xi) * (1.0f / AzimuthalResolution);
        pdf = _sampler->pdf(v, x) * float(AzimuthalResolution * Constant::INV_TWO_PI);

        if(isnan(phi)){
            _sampler->warp(v, xi, x);
        }
    }

