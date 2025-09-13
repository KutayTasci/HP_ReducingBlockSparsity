#include "HPNBM/junk.h"
#include <cstdlib> // For rand()
#include <iostream> // For std::cout
#include <vector>


struct LLNode {
    size_t row;
    size_t col;
    LLNode* next;
};

void generateSparseMAtrix(int rows, int cols, double sparsity, size_t*& ia, size_t*& ja) {
    srand(42);

    // Estimate nnz for allocation (upper bound)
    size_t max_nnz = static_cast<size_t>(rows * cols * sparsity * 1.2) + 1;
    std::vector<size_t> row_counts(rows, 0);
    std::vector<size_t> temp_ja;
    temp_ja.reserve(max_nnz);

    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if ((double)rand() / RAND_MAX <= sparsity) {
                temp_ja.push_back(j);
                row_counts[i]++;
            }
        }
    }

    size_t nnz = temp_ja.size();
    ia = new size_t[rows + 1];
    ja = new size_t[nnz];

    ia[0] = 0;
    for (size_t i = 0; i < rows; ++i) {
        ia[i + 1] = ia[i] + row_counts[i];
    }

    // Fill ja array row-wise
    size_t idx = 0;
    for (size_t i = 0; i < rows; ++i) {
        size_t count = row_counts[i];
        for (size_t k = 0; k < count; ++k) {
            ja[idx++] = temp_ja[idx];
        }
    }
}
