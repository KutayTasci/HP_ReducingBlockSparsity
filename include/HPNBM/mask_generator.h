#ifndef HPNBM_MASK_GENERATOR_H
#define HPNBM_MASK_GENERATOR_H
#include <cstddef>
#include <vector>

void generate_mask(size_t rows, size_t cols, int external_global, int internal_global, 
    int sliding_window, int dilation, int random_per_row, std::vector<size_t> segment_sizes, 
    std::vector<size_t> dilation_sizes, bool apply_causal, size_t*& ia, size_t*& ja);
#endif // HPNBM_MASK_GENERATOR_H