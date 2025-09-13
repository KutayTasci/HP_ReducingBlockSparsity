#ifndef HPNBM_REORDER_H
#define HPNBM_REORDER_H

#include <cstddef>
#include "HPNBM/util.h"


/**
 * @brief Reorders a sparse matrix using the Reverse Cuthill-McKee (RCM) algorithm and converts it to BlockCSR format.
 *
 * This function takes the input sparse matrix in CSR format (represented by ia and ja arrays),
 * applies the RCM algorithm to reduce the matrix bandwidth, and then partitions the reordered matrix
 * into blocks of the specified size, storing the result in a BlockCSR structure.
 *
 * @param ia Pointer to the row pointer array of the input CSR matrix (size n+1).
 * @param ja Pointer to the column indices array of the input CSR matrix.
 * @param n Number of rows (and columns) in the square matrix.
 * @param block_size Size of the blocks for BlockCSR conversion.
 * @param bcsr Reference to a pointer where the resulting BlockCSR matrix will be stored.
 * @param verbose If true, enables verbose output for debugging and progress information.
 */
void reorder_RCM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose=false);
/*
 * @brief Reorders a sparse matrix using the Hypergraph Partititoning Single Border algorithm and converts it to BlockCSR format.

 */
void reorder_HPSB(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose=false);
/*
 * @brief Reorders a sparse matrix using the Hypergraph Partititoning Non-empty block minimization algorithm and converts it to BlockCSR format.

 */
void reorder_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose=false);

/*
 * @brief Reorders a sparse matrix using the combination of Reverse Cuthill-McKee and Hypergraph Partititoning Non-empty block minimization algorithms and converts it to BlockCSR format.

 */
void reorder_RCM_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose=false);
/*
 * @brief Reorders a sparse matrix using the combination of Hypergraph Partititoning Single Border and Hypergraph Partititoning Non-empty block minimization algorithms and converts it to BlockCSR format.

 */
void reorder_HPSB_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose=false);
#endif // HPNBM_REORDER_H