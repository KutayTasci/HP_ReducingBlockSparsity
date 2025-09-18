#ifndef UTIL_H
#define UTIL_H

#include <string>


/**
 * @brief BlockCSR is a data structure representing a block compressed sparse row (CSR) matrix format.
 *
 * This structure is designed to efficiently store and manipulate sparse matrices partitioned into blocks.
 * It contains arrays and metadata for block-wise and row/column-wise indexing, mapping, and permutation.
 *
 * Members:
 * - block_indices: Pointer to array of size (num_block + 1), marking the start of each block in the CSR structure.
 * - ia: Pointer to array representing the row pointer for CSR format within blocks.
 * - ja: Pointer to array representing the column indices for CSR format within blocks.
 * - row_map: Pointer to array mapping global row indices to block-local row indices.
 * - block_id_map: Pointer to array mapping each row to its corresponding block ID.
 * - block_row_ind: Pointer to array containing the row indices of each block.
 * - block_col_ind: Pointer to array containing the column indices of each block.
 * - row_perm: Pointer to array representing a permutation of the rows (useful for reordering).
 * - col_perm: Pointer to array representing a permutation of the columns (useful for reordering).
 * - num_block_rows: Number of block rows in the matrix.
 * - num_block_cols: Total number of blocks (including empty ones).
 * - num_blocks: Number of non-empty blocks in the matrix.
 *
 * This structure is useful for block-based sparse matrix operations, such as block matrix multiplication,
 * and for algorithms that benefit from block partitioning and reordering.
 */
struct BlockCSR {
    size_t* block_indices;    // Size: num_block + 1
    size_t* ia;      
    size_t* ja;      
    size_t* row_map;          // Size: total number of rows
    size_t* block_id_map;
    size_t* block_row_ind;    // Size: num_block_rows + 1
    size_t* block_col_ind;    // Size: num_block_cols + 1
    size_t* row_perm;         // permutation array for rows
    size_t* col_perm;         // permutation array for columns
    size_t num_block_rows;
    size_t num_block_cols;    // total number of blocks (including empty ones)
    size_t num_blocks;        // number of non-empty blocks
};

void testBlockSparsity(size_t* block_row_ptr, size_t* block_col_ind, size_t num_blocks, size_t* ia, size_t* ja, size_t rows, size_t cols);

void ReverseCuthillMcKee(size_t* ia, size_t* ja, size_t n, size_t*& perm);

void reorderCSR(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* row_perm, size_t* col_perm, size_t*& ia_new, size_t*& ja_new);
void reorderCSR_Rows(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* row_perm, size_t*& ia_new, size_t*& ja_new);
void reorderCSR_Cols(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* col_perm, size_t*& ia_new, size_t*& ja_new);

void create_BlockCSR(size_t* block_row_ind, size_t* block_col_ind, size_t num_blocks, size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* row_perm, size_t* column_perm, BlockCSR*& bcsr);
void analyzeBlockCSR(BlockCSR* bcsr);

#endif