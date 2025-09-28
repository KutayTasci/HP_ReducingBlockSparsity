#ifndef HP_H
#define HP_H

#include <cstddef>


void HPSB_RowNet(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t num_parts, size_t*& row_part, size_t*& col_part);
void HPNBM(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t num_parts, size_t block_size, size_t*& row_part, size_t*& block_row_ptr);

void HPNBM_Patoh(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t num_parts, size_t block_size, size_t*& row_part, size_t*& block_row_ptr);

void HP_Rownet(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t num_parts, size_t block_size, size_t*& perm, size_t*& block_row_ptr, size_t*& block_col_ptr);
void HP_TwoConstraint(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t num_parts, size_t block_size, size_t*& row_perm, size_t*& col_perm, size_t*& block_row_ptr, size_t*& block_col_ptr);


#endif