#include "HPNBM/reorder.h"
#include "HPNBM/util.h"
#include "HPNBM/hp.h"
#include "HPNBM/mm.h"


#include <vector>
#include <ctime>
#include <chrono>
#include <iostream>
#include <cstdlib> 
#include <cstring>


std::string write_path = "./tmp/";

void set_write_path(const std::string& path) {
    write_path = path;
}

void writeReorderedMatrix(size_t* ia, size_t* ja, size_t rows, size_t cols, const std::string& filename) {
    std::string out_file = write_path + filename;
    std::cout << "Writing reordered matrix to: " << out_file << std::endl;
    char * cstr = new char[out_file.length() + 1];
    std::strcpy(cstr, out_file.c_str());
    csr2mm(ia, ja, rows, cols, cstr);
    delete[] cstr;
}

void reorder_baseline(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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
    std::cout << "No reordering applied." << std::endl;
    writeReorderedMatrix(ia, ja, rows, cols, "original_matrix.mtx");

    create_BSR(ia, ja, rows, cols, block_size, bsr);

    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details without reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;


}

void reorder_RCM(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR(ia, ja, rows, cols, perm, perm, ia_new, ja_new);
    
    writeReorderedMatrix(ia_new, ja_new, rows, cols, "rcm_reordered.mtx");

    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);

    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after RCM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }
    

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] perm;

    delete[] ia_new;
    delete[] ja_new;

}

void reorder_HPSB(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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
    size_t* row_perm = nullptr;
    printf("Starting HPSB_RowNet...\n");
    std::chrono::system_clock::time_point start_time_hpsb = std::chrono::system_clock::now();
    HPSB_RowNet(ia, ja, rows, cols, number_of_blocks, row_perm, col_perm);

    std::chrono::system_clock::time_point end_time_hpsb = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpsb = end_time_hpsb - start_time_hpsb;
    if (verbose) {
        std::cout << "HPSB_RowNet execution time: " << elapsed_time_hpsb.count() << " seconds." << std::endl;
    }

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    
    reorderCSR(ia, ja, rows, cols, row_perm, col_perm, ia_new, ja_new);
    
    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);
    
    writeReorderedMatrix(ia_new, ja_new, rows, cols, "hpsb_reordered.mtx");

    if (verbose)
    {
        double block_sparsity =  (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPSB reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;

    delete[] col_perm;
    delete[] ia_new;
    delete[] ja_new;

}

void reorder_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR_Rows(ia, ja, rows, cols, row_perm, ia_new, ja_new);
    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);

    writeReorderedMatrix(ia_new, ja_new, rows, cols, "hpnbm_reordered.mtx");
    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] col_perm;
    delete[] row_perm;
    delete[] ia_new;
    delete[] ja_new;
}

void reorder_HPNBM_PaToH(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR_Rows(ia, ja, rows, cols, row_perm, ia_new, ja_new);

    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);

    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPNBM_Patoh reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] col_perm;
    delete[] row_perm;
    delete[] ia_new;
    delete[] ja_new;
}

void reorder_HPRownet(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR_Cols(ia, ja, rows, cols, col_perm, ia_new, ja_new);

    writeReorderedMatrix(ia_new, ja_new, rows, cols, "hprownet_reordered.mtx");

    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);

    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Block sparsity after HP_Rownet: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] col_perm;
    delete[] row_perm;
    delete[] ia_new;
    delete[] ja_new;
    
}

void reorder_RCM_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR_Rows(ia_reordered, ja_reordered, rows, cols, row_perm, ia_new, ja_new);
    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);

    writeReorderedMatrix(ia_new, ja_new, rows, cols, "rcm_hpnbm_reordered.mtx");

    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after RCM+HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] ia_reordered;
    delete[] ja_reordered;
    delete[] row_perm;
    delete[] col_perm;
    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] ia_new;
    delete[] ja_new;

}

void reorder_HPSB_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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
    size_t* row_perm_hpsb = nullptr;
    std::chrono::system_clock::time_point start_time_hpsb = std::chrono::system_clock::now();
    HPSB_RowNet(ia, ja, rows, cols, number_of_blocks, row_perm_hpsb, col_perm_hpsb);


    size_t* ia_reordered = nullptr;
    size_t* ja_reordered = nullptr;

    reorderCSR(ia, ja, rows, cols, row_perm_hpsb, col_perm_hpsb, ia_reordered, ja_reordered);


    size_t* row_perm = nullptr;

    HPNBM(ia_reordered, ja_reordered, rows, cols, number_of_blocks, block_size, row_perm, block_row_ind);
    std::chrono::system_clock::time_point end_time_hpnm = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time_hpnm = end_time_hpnm - start_time_hpsb;

    if (verbose) {
        std::cout << "HPSB+HPNBM execution time: " << elapsed_time_hpnm.count() << " seconds." << std::endl;
    }

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR_Rows(ia_reordered, ja_reordered, rows, cols, row_perm, ia_new, ja_new);
    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);
    writeReorderedMatrix(ia_new, ja_new, rows, cols, "hpsb_hpnbm_reordered.mtx");
    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HPSB+HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] ia_new;
    delete[] ja_new;

    delete[] ia_reordered;
    delete[] ja_reordered;
    delete[] row_perm;

}

void reorder_HPRownet_HPNBM(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR_Rows(ia_reordered, ja_reordered, rows, cols, row_perm, ia_new, ja_new);
    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);
    writeReorderedMatrix(ia_new, ja_new, rows, cols, "hp_hpnbm_reordered.mtx");
    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after HP+HPNBM reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] ia_reordered;
    delete[] ja_reordered;
    delete[] row_perm;
    delete[] col_perm_hp;
    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] ia_new;
    delete[] ja_new;

}

void reorder_TwoConstraint(size_t* ia, size_t* ja, size_t n, size_t block_size, BSR*& bsr, bool verbose) {
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

    size_t* ia_new = nullptr;
    size_t* ja_new = nullptr;
    reorderCSR(ia, ja, rows, cols, row_perm, col_perm, ia_new, ja_new);
    create_BSR(ia_new, ja_new, rows, cols, block_size, bsr);
    writeReorderedMatrix(ia_new, ja_new, rows, cols, "two_constraint_hp_reordered.mtx");
    if (verbose)
    {
        double block_sparsity = (double)(bsr->num_blocks) / (number_of_blocks * number_of_blocks);
        std::cout << "Sparsity details after Two-Constraint HP reordering: " << block_sparsity * 100 << "%" << std::endl;
    }

    delete[] block_row_ind;
    delete[] block_col_ind;
    delete[] col_perm;
    delete[] row_perm;
    delete[] ia_new;
    delete[] ja_new;

}
