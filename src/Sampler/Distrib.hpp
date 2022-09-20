#pragma once
#include "../Common/math.hpp"
#include "vector"
struct Distribution1D{

    Distribution1D(const Float *f, int n) : func(f, f + n), cdf(n + 1) {
        // Compute integral of step function at $x_i$
        cdf[0] = 0;
        for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[n];
        if (funcInt == 0) {
            for (int i = 1; i < n + 1; ++i) cdf[i] = Float(i) / Float(n);
        } else {
            for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
        }
    }

    int SampleDiscrete(Float sample,Float * pdf) const {
        int offset = FindInterval((int)cdf.size(),
                                           [&](int index) { return cdf[index] <= sample; });
        if(pdf) *pdf= DiscretePDF(offset);
        return offset;
    }

    Float SampleContinuous(Float u, Float *pdf, int *off = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval((int)cdf.size(),
                                  [&](int index) { return cdf[index] <= u; });
        if (off) *off = offset;
        // Compute offset along CDF segment
        Float du = u - cdf[offset];
        if ((cdf[offset + 1] - cdf[offset]) > 0) {
            du /= (cdf[offset + 1] - cdf[offset]);
        }

        // Compute PDF for sampled offset
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / Count();
    }

    Float DiscretePDF(int index) const {
        if(index>= func.size()){
            return 0;
        }
        return func[index] / (funcInt * Count());
    }

    void warp(Float & sample,int  & index){
        index = FindInterval((int)cdf.size(),
                             [&](int index) { return cdf[index] <= sample; });
        sample = (sample-cdf[index])/func[index];
    }

    int Count() const {
        return int(func.size());
    }
    std::vector<Float> func, cdf;
    Float funcInt;
};


struct Distribution2D{
public:
    // Distribution2D Public Methods
    Distribution2D(const Float *data, int nu, int nv);
    vec2 SampleContinuous(const vec2 &u, Float *pdf) const {
        Float pdfs[2];
        int v;
        Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
        Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
        if(pdf)
        *pdf = pdfs[0] * pdfs[1];
        return vec2(d0, d1);
    }
    Float Pdf(const vec2 &p) const {
        int iu = clamp(int(p[0] * pConditionalV[0]->Count()), 0,
                       pConditionalV[0]->Count() - 1);
        int iv = clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
        return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
    }

    void  warp(vec2 &uv, int &row, int &column) const;

private:
    // Distribution2D Private Data
    std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
    std::unique_ptr<Distribution1D> pMarginal;
};