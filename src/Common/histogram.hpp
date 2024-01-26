#pragma once

#include <vector>

class Histogram {
public:
    Histogram(const std::vector<float>& data, size_t num_bins);

    float level(float count_percentage) const;

    std::vector<size_t> counts;
    float               bin_size;
    size_t              data_size;
};