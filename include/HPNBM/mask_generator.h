#ifndef HPNBM_MASK_GENERATOR_H
#define HPNBM_MASK_GENERATOR_H
#include <cstddef>
#include <vector>

void generate_mask(size_t rows, size_t cols, int external_global, int internal_global, 
    int sliding_window, int dilation, int random_per_row, std::vector<size_t> segment_sizes, 
    std::vector<size_t> dilation_sizes, bool apply_causal, size_t*& ia, size_t*& ja);

void generate_bigbird_mask(
    size_t rows,
    size_t cols,
    int window_size,
    int num_global_tokens,
    int random_per_row,
    bool apply_causal, 
    size_t*& ia,
    size_t*& ja
);

void generate_longformer_mask(
    size_t rows,
    size_t cols,
    int window_size,          // local window
    int dilation,              // dilation for sliding window
    int num_internal_global_tokens,    // global tokens
    bool apply_causal,        // optionally causal
    size_t*& ia,
    size_t*& ja
);

void generate_longnet_mask(
    size_t rows,
    size_t cols,
    int window_size,                                      // local window (LongNet still has this)
    const std::vector<size_t>& segment_sizes,             // user-provided segments
    const std::vector<size_t>& dilation_sizes,            // user-provided dilations
    bool apply_causal,                                    // optionally causal
    size_t*& ia,
    size_t*& ja
);
#endif // HPNBM_MASK_GENERATOR_H