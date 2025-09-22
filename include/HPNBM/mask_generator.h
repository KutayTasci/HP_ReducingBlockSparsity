#ifndef HPNBM_MASK_GENERATOR_H
#define HPNBM_MASK_GENERATOR_H
#include <cstddef>

void generate_mask(size_t rows, size_t cols, int external_global, int internal_global, int sliding_window, int dilation, int random_per_row, bool apply_causal, size_t*& ia, size_t*& ja);
#endif // HPNBM_MASK_GENERATOR_H