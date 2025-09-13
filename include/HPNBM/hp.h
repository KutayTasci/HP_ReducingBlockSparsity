#ifndef HP_H
#define HP_H

#include <cstddef>


void HPSB_RowNet(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t num_parts, size_t*& col_part);
void HPNBM(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t num_parts, size_t block_size, size_t*& row_part, size_t*& block_row_ptr);

#endif