#include "HPNBM/reorder.h"
#include "HPNBM/util.h"
#include "HPNBM/hp.h"


#include <vector>
#include <ctime>
#include <chrono>
#include <iostream>
#include <cstdlib> 


void reorder_RCM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];


    size_t rows = n;
    size_t cols = n;

    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }


    
    size_t* perm = nullptr;
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    ReverseCuthillMcKee(ia, ja, rows, perm);
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_rcm = end - start;
    if (verbose) {
        std::cout << "ReverseCuthillMcKee execution time: " << elapsed_time_rcm.count() << " seconds." << std::endl;
    }


    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols, perm, perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after RCM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }
    

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] perm;

}

void reorder_HPSB(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;
    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* col_perm = nullptr;
    std::chrono::system_clock::time_point start_time_hpsb = std::chrono::system_clock::now();
    HPSB_RowNet(ia, ja, rows, cols, number_of_blocks, col_perm);

    std::chrono::system_clock::time_point end_time_hpsb = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpsb = end_time_hpsb - start_time_hpsb;
    if (verbose) {
        std::cout << "HPSB_RowNet execution time: " << elapsed_time_hpsb.count() << " seconds." << std::endl;
    }

    size_t* row_perm = new size_t[rows];
    for (size_t i = 0; i < rows; i++) {
        row_perm[i] = i;
    }

    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity =  (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPSB reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;

    delete[] col_perm;
    delete[] row_perm;

}

void reorder_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;
    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* row_perm = nullptr;
    size_t* col_perm = nullptr;

    std::chrono::system_clock::time_point start_time_hpnm = std::chrono::system_clock::now();
    HPNBM(ia, ja, rows, cols, number_of_blocks, block_size, row_perm, block_row_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hpnm;

    if (verbose) {
        std::cout << "HPSB_RowNet execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    col_perm = new size_t[rows];
    for (size_t i = 0; i < rows; i++) {
        col_perm[i] = i;
    }

    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] col_perm;
    delete[] row_perm;
}

void reorder_HPNBM_PaToH(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;
    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* row_perm = nullptr;
    size_t* col_perm = nullptr;

    std::chrono::system_clock::time_point start_time_hpnm = std::chrono::system_clock::now();
    HPNBM_Patoh(ia, ja, rows, cols, number_of_blocks, block_size, row_perm, block_row_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hpnm;

    if (verbose) {
        std::cout << "HPSB_RowNet_Patoh execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    col_perm = new size_t[rows];
    for (size_t i = 0; i < rows; i++) {
        col_perm[i] = i;
    }

    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPNBM_Patoh reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] col_perm;
    delete[] row_perm;
}

void reorder_HPRownet(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;

    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* row_perm = nullptr;
    size_t* col_perm = nullptr;

    std::chrono::system_clock::time_point start_time_hpnm = std::chrono::system_clock::now();
    HP_Rownet(ia, ja, rows, cols, number_of_blocks, block_size, col_perm, block_row_ind, block_col_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hpnm;

    if (verbose) {
        std::cout << "HP_Rownet execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    row_perm = new size_t[cols];
    for (size_t i = 0; i < cols; i++) {
        row_perm[i] = i;
    }

    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Block sparsity after HP_Rownet: " << block_sparsity * 100 << "%" << std::endl;
    }
    
}

void reorder_RCM_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;

    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* perm = nullptr;
    std::chrono::system_clock::time_point start_time_hpsb = std::chrono::system_clock::now();
    ReverseCuthillMcKee(ia, ja, rows, perm);


    size_t* ia_reordered = nullptr;
    size_t* ja_reordered = nullptr;
    reorderCSR(ia, ja, rows, cols, perm, perm, ia_reordered, ja_reordered);

    size_t* row_perm = perm;

    size_t* col_perm = nullptr;

    HPNBM(ia_reordered, ja_reordered, rows, cols, number_of_blocks, block_size, row_perm, block_row_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hpsb;

    if (verbose) {
        std::cout << "RCM+HPNBM execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    col_perm = new size_t[cols];
    for (size_t i = 0; i < cols; i++) {
        col_perm[i] = i;
    }
    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia_reordered, ja_reordered, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after RCM+HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] ia_reordered;
    delete[] ja_reordered;
    delete[] row_perm;
    delete[] col_perm;

}

void reorder_HPSB_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;

    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* col_perm_hpsb = nullptr;
    std::chrono::system_clock::time_point start_time_hpsb = std::chrono::system_clock::now();
    HPSB_RowNet(ia, ja, rows, cols, number_of_blocks, col_perm_hpsb);


    size_t* ia_reordered = nullptr;
    size_t* ja_reordered = nullptr;

    reorderCSR_Cols(ia, ja, rows, cols, col_perm_hpsb, ia_reordered, ja_reordered);


    size_t* row_perm = nullptr;

    HPNBM(ia_reordered, ja_reordered, rows, cols, number_of_blocks, block_size, row_perm, block_row_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hpsb;

    if (verbose) {
        std::cout << "HPSB+HPNBM execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    size_t* col_perm = new size_t[cols];
    for (size_t i = 0; i < cols; i++) {
        col_perm[i] = i;
    }
    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia_reordered, ja_reordered, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPSB+HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    

    delete[] ia_reordered;
    delete[] ja_reordered;
    delete[] row_perm;
    delete[] col_perm;

}

void reorder_HPRownet_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;

    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* col_perm_hp = nullptr;
    std::chrono::system_clock::time_point start_time_hp = std::chrono::system_clock::now();
    HP_Rownet(ia, ja, rows, cols, number_of_blocks, block_size, col_perm_hp, block_row_ind, block_col_ind);

    size_t* ia_reordered = nullptr;
    size_t* ja_reordered = nullptr;

    reorderCSR_Cols(ia, ja, rows, cols, col_perm_hp, ia_reordered, ja_reordered);

    size_t* row_perm = nullptr;

    HPNBM(ia_reordered, ja_reordered, rows, cols, number_of_blocks, block_size, row_perm, block_row_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hp;

    if (verbose) {
        std::cout << "HP+HPNBM execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    size_t* col_perm = new size_t[cols];
    for (size_t i = 0; i < cols; i++) {
        col_perm[i] = i;
    }
    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia_reordered, ja_reordered, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HP+HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] ia_reordered;
    delete[] ja_reordered;
    delete[] row_perm;
    delete[] col_perm;
}

void reorder_TwoConstraint(size_t* ia, size_t* ja, size_t n, size_t block_size, BlockCSR*& bcsr, bool verbose) {
    size_t number_of_blocks = n / block_size;
    size_t* block_row_ind = new size_t[number_of_blocks + 1];
    size_t* block_col_ind = new size_t[number_of_blocks + 1];

    size_t rows = n;
    size_t cols = n;

    block_row_ind[0] = 0;
    for (size_t i = 1; i <= number_of_blocks; ++i) {
        block_row_ind[i] = i * block_size;
    }

    block_col_ind[0] = 0;
    for (size_t j = 1; j <= number_of_blocks; ++j) {
        block_col_ind[j] = j * block_size;
    }

    if (verbose) {
        std::cout << "---------------------------------" << std::endl;
        testBlockSparsity(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols);
    }

    size_t* row_perm = nullptr;
    size_t* col_perm = nullptr;

    std::chrono::system_clock::time_point start_time_hpnm = std::chrono::system_clock::now();
    HP_TwoConstraint(ia, ja, rows, cols, number_of_blocks, block_size, row_perm, col_perm, block_row_ind, block_col_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hpnm;

    if (verbose) {
        std::cout << "Two-Constraint HP execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    create_BlockCSR(block_row_ind, block_col_ind, number_of_blocks, ia, ja, rows, cols, row_perm, col_perm, bcsr);

    if (verbose)
    {
        double block_sparsity = (double)(bcsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after Two-Constraint HP reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

}
