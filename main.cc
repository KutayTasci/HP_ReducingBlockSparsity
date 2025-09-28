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
#include "HPNBM/mm.h"



extern char* write_path;

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

    BSR* bsr = nullptr;

    std::cout << "========================" << std::endl;
    std::cout << "Experimental results for Baseline" << std::endl;
    std::cout << "========================" << std::endl;
    reorder_baseline(ia, ja, rows, block_size, bsr, false);
    analyzeBSR(bsr);

    freeBSR(bsr);
    

    if (RCM) {
        BSR* bsr_rcm = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for RCM" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_RCM(ia, ja, rows, block_size, bsr_rcm, false);
        analyzeBSR(bsr_rcm);

        freeBSR(bsr_rcm);


    }

    if (HPRownet) {
        BSR* bsr_hpr = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Rownet_HyperGraph_Partitioning" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPRownet(ia, ja, rows, block_size, bsr_hpr, false);
        analyzeBSR(bsr_hpr);
        freeBSR(bsr_hpr);
    }

    if (HPSB) {
        BSR* bsr_hpsb = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Single_Border" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPSB(ia, ja, rows, block_size, bsr_hpsb, false);
        analyzeBSR(bsr_hpsb);
        freeBSR(bsr_hpsb);
    }

    if (HPNBM) {
        BSR* bsr_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPNBM(ia, ja, rows, block_size, bsr_hpn, false);
        analyzeBSR(bsr_hpn);
        freeBSR(bsr_hpn);
    }

    if (TwoConstraint) {
        BSR* bsr_two = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Two_Constraint" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_TwoConstraint(ia, ja, rows, block_size, bsr_two, false);
        analyzeBSR(bsr_two);
        freeBSR(bsr_two);
    }

    if (RCM_HPNBM) {
        BSR* bsr_rcm_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for RCM + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_RCM_HPNBM(ia, ja, rows, block_size, bsr_rcm_hpn, false);
        analyzeBSR(bsr_rcm_hpn);
        freeBSR(bsr_rcm_hpn);
    }

    if (HPSB_HPNBM) {
        BSR* bsr_hpsb_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Single_Border + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPSB_HPNBM(ia, ja, rows, block_size, bsr_hpsb_hpn, false);
        analyzeBSR(bsr_hpsb_hpn);
        freeBSR(bsr_hpsb_hpn);
    }

    if (HPRownet_HPNBM) {
        BSR* bsr_hpr_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Rownet_HyperGraph_Partitioning + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPRownet_HPNBM(ia, ja, rows, block_size, bsr_hpr_hpn, false);
        analyzeBSR(bsr_hpr_hpn);
        freeBSR(bsr_hpr_hpn);
    }
}

void run_experiments() {  

    // Save the original cout buffer
    std::streambuf* original_cout_buffer = std::cout.rdbuf();

    // Create an output file stream
    std::ofstream file_out("output.txt");

    // Redirect cout to the file stream
    std::cout.rdbuf(file_out.rdbuf());
 
    bool RCM = false;
    bool HPRownet = false;
    bool HPSB = true;
    bool HPNBM = false;
    bool TwoConstraint = false;
    bool RCM_HPNBM = false;
    bool HPSB_HPNBM = true;
    bool HPRownet_HPNBM = false;

    int num_experiments = 1;

    size_t list_of_number_of_blocks[] = {1024}; //1024, 4096, 16384, 65536, 2048, 4096};
    for (size_t number_of_blocks : list_of_number_of_blocks) {
        size_t block_size = 16;
        size_t rows = number_of_blocks * block_size;
        size_t cols = number_of_blocks * block_size;
        std::cout << "========================" << std::endl;
        std::cout << "Running experiments for matrix size: " << rows << " x " << cols << std::endl;
        std::cout << "========================" << std::endl;
        for (int exp = 0; exp < num_experiments; ++exp) {
            // write as char array
            write_path = "./tmp/";
            // create directory if it does not exist
            std::string command = "mkdir -p " + std::string(write_path);
            system(command.c_str());
            printf("Experiment %d of %d\n", exp + 1, num_experiments);
            std::cout << "========================" << std::endl;
            std::cout << "Experiment " << (exp + 1) << " of " << num_experiments << std::endl;
            std::cout << "Matrix size: " << rows << " x " << cols << std::endl;
            std::cout << "Block size: " << block_size << " x " << block_size << std::endl;
            size_t* ia = nullptr;
            size_t* ja = nullptr;
            int seg_len = std::log2(rows) - 2;
            std::vector<size_t> segment_sizes(seg_len);
            std::vector<size_t> dilation_sizes(seg_len);

            for (int i = 0; i < seg_len; ++i) {
                segment_sizes[i] = std::pow(2, i+2);
                dilation_sizes[i] = std::pow(2, i);
            }
            segment_sizes = {};
            dilation_sizes = {};
            double sparsity = 0.001;
            int external_global = 0; //std::log2(rows);
            int internal_global = 0; //std::log2(rows);
            int sliding_window = 4 ; //std::log2(rows);
            int sliding_dilation = 0;
            int random_per_row = 3; //sparsity * cols ; //std::log2(rows);
            bool causal = true;
            generate_mask(rows, cols, external_global, internal_global, sliding_window, sliding_dilation, random_per_row, segment_sizes, dilation_sizes, causal, ia, ja);

            std::cout << "========================" << std::endl;
            
            run_experiments(RCM, HPRownet, HPSB, HPNBM, TwoConstraint, RCM_HPNBM, HPSB_HPNBM, HPRownet_HPNBM, ia, ja, rows, block_size);
            //clean up
            free(ia);
            free(ja);
            segment_sizes.clear();
            dilation_sizes.clear();
        }
    }
}

void run_with_mm(char* filename) {
    mmdata* mm = (struct mmdata *) calloc(1, sizeof(struct mmdata));
    std::cout << "========================" << std::endl;
    if (initialize_mm(filename, mm) != 0) {
        std::cerr << "Error initializing matrix from file: " << filename << std::endl;
        return;
    }

    
    size_t* ia = nullptr;
    size_t* ja = nullptr;
    size_t rows, cols;
    mm2csr(mm, ia, ja, rows, cols);

    //printCSR(ia, ja, rows);

    size_t block_size = 16; // Example block size
    run_experiments(true, true, true, true, true, true, true, true, ia, ja, rows, block_size);

    freemm(mm);
    free(ia);
    free(ja);
}

void reorder_and_save() {

}

int main(int argc, char* argv[]) {

    // Run experiments
    run_experiments();

    //run_with_mm("../sparse_matrix_groot.mtx");

    return 0; // Indicate successful execution
}