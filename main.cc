#include <iostream>   // For input/output
#include <string>     // For using std::string
#include <vector>     // For using std::vector (if needed)
#include <ctime>      // For measuring execution time
#include <cmath>
#include <fstream>
#include <streambuf> // Required for std::streambuf

#include "HPNBM/util.h"
#include "HPNBM/reorder.h"
#include "HPNBM/mask_generator.h"



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

void run_experiments(bool RCM, bool HPRownet, bool HPSB, bool HPNBM, bool TwoConstraint, 
    bool RCM_HPNBM, bool HPSB_HPNBM, bool HPRownet_HPNBM, size_t* ia, size_t* ja, size_t rows, size_t block_size) {
    
    BlockCSR* bcsr = nullptr;
    
    std::cout << "========================" << std::endl;
    std::cout << "Experimental results for Baseline" << std::endl;
    std::cout << "========================" << std::endl;
    reorder_baseline(ia, ja, rows, block_size, bcsr, false);
    analyzeBlockCSR(bcsr);
    bcsr = nullptr;

    if (RCM) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for RCM" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_RCM(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }

    if (HPRownet) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Rownet_HyperGraph_Partitioning" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPRownet(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }

    if (HPSB) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Single_Border" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPSB(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }

    if (HPNBM) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPNBM(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }

    if (TwoConstraint) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Two_Constraint" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_TwoConstraint(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }

    if (RCM_HPNBM) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for RCM + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_RCM_HPNBM(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }

    if (HPSB_HPNBM) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Single_Border + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPSB_HPNBM(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }

    if (HPRownet_HPNBM) {
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Rownet_HyperGraph_Partitioning + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPRownet_HPNBM(ia, ja, rows, block_size, bcsr, false);
        analyzeBlockCSR(bcsr);
        bcsr = nullptr;
    }
}

int main(int argc, char* argv[]) {
    
    
    


    //printf("Generating sparse matrix with %zu rows, %zu cols, and %.5f sparsity.\n", rows, cols, sparsity);
    //generateSparseMAtrix(rows, cols, sparsity, ia, ja);
    

    
    //printCSR(ia, ja, rows);
    //printf("Sparse matrix generated.\n");
    
    // Save the original cout buffer
    std::streambuf* original_cout_buffer = std::cout.rdbuf();

    // Create an output file stream
    std::ofstream file_out("output.txt");

    // Redirect cout to the file stream
    std::cout.rdbuf(file_out.rdbuf());
 
    bool RCM = true;
    bool HPRownet = true;
    bool HPSB = true;
    bool HPNBM = true;
    bool TwoConstraint = true;
    bool RCM_HPNBM = true;
    bool HPSB_HPNBM = true;
    bool HPRownet_HPNBM = true;

    int num_experiments = 2;

    size_t list_of_number_of_blocks[] = {256, 512, 1024};
    for (size_t number_of_blocks : list_of_number_of_blocks) {
        size_t block_size = 16;
        size_t rows = number_of_blocks * block_size;
        size_t cols = number_of_blocks * block_size;
        std::cout << "========================" << std::endl;
        std::cout << "Running experiments for matrix size: " << rows << " x " << cols << std::endl;
        std::cout << "========================" << std::endl;
        for (int exp = 0; exp < num_experiments; ++exp) {
            std::cout << "========================" << std::endl;
            std::cout << "Experiment " << (exp + 1) << " of " << num_experiments << std::endl;
            std::cout << "Matrix size: " << rows << " x " << cols << std::endl;
            std::cout << "Block size: " << block_size << " x " << block_size << std::endl;
            size_t* ia = nullptr;
            size_t* ja = nullptr;
            int seg_len = std::log2(rows) - 2;
            size_t* segment_sizes = new size_t[seg_len];
            size_t* dilation_sizes = new size_t[seg_len];
            printf("Segment length: %d\n", seg_len);
            for (int i = 0; i < seg_len; ++i) {
                segment_sizes[i] = std::pow(2, i+2);
                dilation_sizes[i] = std::pow(2, i);
            }
            int external_global = std::log2(rows);
            int internal_global = std::log2(rows);
            int sliding_window = std::log2(rows);
            int sliding_dilation = 0;
            int random_per_row = std::log2(rows);
            bool causal = true;
            generate_mask(rows, cols, external_global, internal_global, sliding_window, sliding_dilation, random_per_row, segment_sizes, dilation_sizes, causal, ia, ja);
            std::cout << "========================" << std::endl;
            
            run_experiments(RCM, HPRownet, HPSB, HPNBM, TwoConstraint, RCM_HPNBM, HPSB_HPNBM, HPRownet_HPNBM, ia, ja, rows, block_size);
        }
    }
    return 0; // Indicate successful execution
}