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



extern std::string write_path;

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

void run_experiments(int no_ofblocks, double sparsity, int seg_start=2) {  
    // ====== User Configurable Parameters ======
    // Set the experiment parameters here
    std::string folder_for_current_run = "/scratch/pioneer/users/kxt437/Dilated2D_FSeg" + std::to_string((int) std::pow(2, seg_start)) + "_N" + std::to_string(no_ofblocks) + "_Causal"+ "/";
    // std::to_string((int) std::pow(2, seg_start))  std::to_string(sparsity)
    // Bigbird_FSeg


    size_t block_size = 16;
    bool causal = true;
    int num_experiments = 10;
    std::vector<size_t> list_of_number_of_blocks = {(size_t)no_ofblocks}; // e.g., {1024, 2048, 4096, 8192, 16384, 32768, 65536}

    int big_bird_blocks = ((no_ofblocks * block_size) * sparsity) / 8;
    // Mask generation parameters (set to 0 or fill as needed)
    int external_global = 0;//2 * big_bird_blocks; // e.g., 2 * big_bird_blocks
    int internal_global = 0;
    int sliding_window = 0;//3 * big_bird_blocks; // e.g., 3 * big_bird_blocks
    int sliding_dilation = 0;
    int random_per_row = 0; //3 * big_bird_blocks; // will be set below
    std::vector<size_t> segment_sizes = {}; // e.g., {4, 8, 16}
    std::vector<size_t> dilation_sizes = {}; // e.g., {1, 2, 4}

    for (int seg = seg_start; seg < static_cast<int>(std::log2(block_size * no_ofblocks)); ++seg) {
        segment_sizes.push_back(static_cast<size_t>(std::pow(2, seg)));
        dilation_sizes.push_back(static_cast<size_t>(std::pow(2, seg - seg_start)));
    }
    //segment_sizes = {};
    //dilation_sizes = {};
    
    // Which experiments to run
    bool RCM = true;
    bool HPRownet = true;
    bool HPSB = true;
    bool HPNBM = false;
    bool TwoConstraint = false;
    bool RCM_HPNBM = false;
    bool HPSB_HPNBM = true;
    bool HPRownet_HPNBM = true;
    // ==========================================

    // Save the original cout buffer
    std::streambuf* original_cout_buffer = std::cout.rdbuf();

    // Create output folder
    std::string command = "mkdir -p " + folder_for_current_run;
    system(command.c_str());

    for (size_t number_of_blocks : list_of_number_of_blocks) {
        std::string output_filename = folder_for_current_run + "/" + std::to_string(number_of_blocks) + "_blocks.txt";
        std::ofstream file_out(output_filename);
        std::cout.rdbuf(file_out.rdbuf());

        size_t rows = number_of_blocks * block_size;
        size_t cols = number_of_blocks * block_size;
        std::cout << "========================" << std::endl;
        std::cout << "Running experiments for matrix size: " << rows << " x " << cols << std::endl;
        std::cout << "========================" << std::endl;

        for (int exp = 0; exp < num_experiments; ++exp) {
            std::string folder_for_current_experiment = folder_for_current_run + "/Exp" + std::to_string(exp) + "/";
            write_path = folder_for_current_experiment;
            std::string command = "mkdir -p " + write_path;
            system(command.c_str());
            set_write_path(write_path);

            printf("Experiment %d of %d\n", exp + 1, num_experiments);
            std::cout << "========================" << std::endl;
            std::cout << "Experiment " << (exp + 1) << " of " << num_experiments << std::endl;
            std::cout << "Matrix size: " << rows << " x " << cols << std::endl;
            std::cout << "Block size: " << block_size << " x " << block_size << std::endl;

            size_t* ia = nullptr;
            size_t* ja = nullptr;

            // Set random_per_row based on sparsity and cols
            //random_per_row = static_cast<int>(sparsity * cols);

            // Generate mask
            generate_mask(rows, cols, external_global, internal_global, sliding_window, sliding_dilation, random_per_row, segment_sizes, dilation_sizes, causal, ia, ja);

            std::cout << "========================" << std::endl;
            run_experiments(RCM, HPRownet, HPSB, HPNBM, TwoConstraint, RCM_HPNBM, HPSB_HPNBM, HPRownet_HPNBM, ia, ja, rows, block_size);

            free(ia);
            free(ja);
        }
    }

    // Restore cout buffer
    std::cout.rdbuf(original_cout_buffer);
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
    int no_ofblocks = argc > 1 ? std::stoi(argv[1]) : 4096;
    double sparsity = argc > 2 ? std::stod(argv[2]) : 0.1;
    int seg_start = argc > 3 ? std::stoi(argv[3]) : 2;
    run_experiments(no_ofblocks, sparsity, seg_start);

    //run_with_mm("../sparse_matrix_groot.mtx");

    return 0; // Indicate successful execution
}