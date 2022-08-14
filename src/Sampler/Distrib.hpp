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
        int idx;
        for(idx=1;cdf[idx]<=Count() && cdf[idx]<sample ;idx++) ;
        idx=idx-1;
        assert(idx<3);
        *pdf= DiscretePDF(idx);
        return idx;
    }

    Float DiscretePDF(int index) const {

        return func[index] / (funcInt * Count());
    }

    int Count() const {
        return int(func.size());
    }
    std::vector<Float> func, cdf;

    Float funcInt;
};