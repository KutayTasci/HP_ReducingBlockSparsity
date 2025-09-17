#include <iostream>   // For input/output
#include <string>     // For using std::string
#include <vector>     // For using std::vector (if needed)
#include <ctime>      // For measuring execution time


#include "HPNBM/util.h"
#include "HPNBM/junk.h"
#include "HPNBM/hp.h"
#include "HPNBM/reorder.h"



void printCSR(size_t* ia, size_t* ja, size_t rows) {
    // Output the generated CSR arrays 

    std::cout << "IA: ";
    for (size_t i = 0; i <= rows; ++i) {
        std::cout << ia[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "JA: ";
    for (size_t i = 0; i < ia[rows]; ++i) {
        std::cout << ja[i] << " ";
    }
    std::cout << std::endl;
   
}

int main(int argc, char* argv[]) {
    

    size_t number_of_blocks = 128;
    size_t block_size = 16;
    size_t rows = number_of_blocks * block_size;
    size_t cols = number_of_blocks * block_size;
    double sparsity = 0.001; 
    size_t* ia = nullptr;
    size_t* ja = nullptr;

    printf("Generating sparse matrix with %zu rows, %zu cols, and %.5f sparsity.\n", rows, cols, sparsity);
    generateSparseMAtrix(rows, cols, sparsity, ia, ja);
    //printCSR(ia, ja, rows);
    printf("Sparse matrix generated.\n");


    BlockCSR* bcsr = nullptr;
    reorder_RCM(ia, ja, rows, block_size, bcsr, true);

    bcsr = nullptr;
    reorder_HPSB(ia, ja, rows, block_size, bcsr, true);

    bcsr = nullptr;
    reorder_HPNBM(ia, ja, rows, block_size, bcsr, true);

    bcsr = nullptr;
    reorder_HPNBM_PaToH(ia, ja, rows, block_size, bcsr, true);


    bcsr = nullptr;
    reorder_HPRownet(ia, ja, rows, block_size, bcsr, true);

    bcsr = nullptr;
    reorder_RCM_HPNBM(ia, ja, rows, block_size, bcsr, true);

    bcsr = nullptr;
    reorder_HPSB_HPNBM(ia, ja, rows, block_size, bcsr, true);


    bcsr = nullptr;
    reorder_HPRownet_HPNBM(ia, ja, rows, block_size, bcsr, true);

    return 0; // Indicate successful execution
}