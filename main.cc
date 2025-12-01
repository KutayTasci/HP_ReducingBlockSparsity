#include <iostream>   // For input/output
#include <string>     // For using std::string
#include <vector>     // For using std::vector (if needed)
#include <ctime>      // For measuring execution time
#include <cmath>
#include <fstream>
#include <streambuf> // Required for std::streambuf
#include <cstring>   // Required for std::strcpy
#include <optional>



#include "HPNBM/util.h"
#include "HPNBM/reorder.h"
#include "HPNBM/mask_generator.h"
#include "HPNBM/mm.h"

//
// Supported mask patterns
//
enum class MaskPattern {
    BigBird,
    LongNet,
    LongFormer
};

//
// Global configuration
//
struct Config {
    int no_blocks = -1;              // required
    int block_size = -1;             // required
    MaskPattern pattern;             // required
    double sparsity = -1.0;          // required for bigbird/longformer
    int segment_start = -1;          // required for longnet
    bool apply_causal = false;       // optional default = false
    std::optional<std::string> save_path = std::nullopt; // optional
    std::optional<std::string> matrix_file = std::nullopt; // optional matrix file
};

static std::string pattern_to_string(MaskPattern p) {
    switch (p) {
        case MaskPattern::BigBird:   return "bigbird";
        case MaskPattern::LongNet:   return "longnet";
        case MaskPattern::LongFormer:return "longformer";
    }
    return "unknown";
}

extern std::string write_path;


// Function to set it
inline void set_global_write_path(const std::string& path) {
    write_path = path;
}

void run_experiments(bool RCM, bool HPRownet, bool HPSB, bool HPNBM, bool TwoConstraint, 
    bool RCM_HPNBM, bool HPSB_HPNBM, bool HPRownet_HPNBM, size_t* ia, size_t* ja, size_t rows, size_t block_size) {

    BSR* bsr = nullptr;

    std::cout << "========================" << std::endl;
    std::cout << "Experimental results for Baseline" << std::endl;
    std::cout << "========================" << std::endl;
    reorder_baseline(ia, ja, rows, block_size, bsr, false);
    //analyzeBSR(bsr);

    //freeBSR(bsr);
    

    if (RCM) {
        BSR* bsr_rcm = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for RCM" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_RCM(ia, ja, rows, block_size, bsr_rcm, false);
        //analyzeBSR(bsr_rcm);

        //freeBSR(bsr_rcm);


    }

    if (HPRownet) {
        BSR* bsr_hpr = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Rownet_HyperGraph_Partitioning" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPRownet(ia, ja, rows, block_size, bsr_hpr, false);
        //analyzeBSR(bsr_hpr);
        //freeBSR(bsr_hpr);
    }

    if (HPSB) {
        BSR* bsr_hpsb = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Single_Border" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPSB(ia, ja, rows, block_size, bsr_hpsb, false);
        //analyzeBSR(bsr_hpsb);
        //freeBSR(bsr_hpsb);
    }

    if (HPNBM) {
        BSR* bsr_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPNBM(ia, ja, rows, block_size, bsr_hpn, false);
        //analyzeBSR(bsr_hpn);
        //freeBSR(bsr_hpn);
    }

    if (TwoConstraint) {
        BSR* bsr_two = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Two_Constraint" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_TwoConstraint(ia, ja, rows, block_size, bsr_two, false);
        //analyzeBSR(bsr_two);
        //freeBSR(bsr_two);
    }

    if (RCM_HPNBM) {
        BSR* bsr_rcm_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for RCM + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_RCM_HPNBM(ia, ja, rows, block_size, bsr_rcm_hpn, false);
        //analyzeBSR(bsr_rcm_hpn);
        //freeBSR(bsr_rcm_hpn);
    }

    if (HPSB_HPNBM) {
        BSR* bsr_hpsb_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Hypergraph_Partitioning_Single_Border + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPSB_HPNBM(ia, ja, rows, block_size, bsr_hpsb_hpn, false);
        //analyzeBSR(bsr_hpsb_hpn);
        //freeBSR(bsr_hpsb_hpn);
    }

    if (HPRownet_HPNBM) {
        BSR* bsr_hpr_hpn = nullptr;
        std::cout << "========================" << std::endl;
        std::cout << "Experimental results for Rownet_HyperGraph_Partitioning + Hypergraph_Partitioning_Nonzero_Block_Minimization" << std::endl;
        std::cout << "========================" << std::endl;
        reorder_HPRownet_HPNBM(ia, ja, rows, block_size, bsr_hpr_hpn, false);
        //analyzeBSR(bsr_hpr_hpn);
        //freeBSR(bsr_hpr_hpn);
    }
}

void run_experiments_mask(const Config& cfg) {
    // Determine pattern string
    std::string pattern_str;
    switch(cfg.pattern) {
        case MaskPattern::BigBird:    pattern_str = "BigBird"; break;
        case MaskPattern::LongNet:    pattern_str = "LongNet"; break;
        case MaskPattern::LongFormer: pattern_str = "LongFormer"; break;
    }

    // Determine folder: use save_path if provided, else temp folder
    std::string folder;
    if (cfg.save_path) {
        folder = cfg.save_path.value();
    } else {
        folder = "tmp/experiments";
    }

    if (pattern_str == "BigBird" || pattern_str == "LongFormer") {
        folder += "/" + pattern_str + "_S" + std::to_string(cfg.sparsity) + "_N" + std::to_string(cfg.no_blocks) + "_B" + std::to_string(cfg.block_size);
    } else if (pattern_str == "LongNet") {
        folder += "/" + pattern_str + "_G" + std::to_string(cfg.segment_start) + "_N" + std::to_string(cfg.no_blocks) + "_B" + std::to_string(cfg.block_size);
    } else {
        std::cout << "Unknown pattern!" << std::endl;
        return;
    }
    if (cfg.apply_causal) folder += "_Causal";
    set_global_write_path(folder);
    
    // Create folder
    std::string command = "mkdir -p " + folder;
    int rs = system(command.c_str());
    if (rs != 0) {
        std::cerr << "Error creating directory: " << folder << std::endl;
        return;
    }

    size_t rows = cfg.no_blocks * cfg.block_size;
    size_t cols = cfg.no_blocks * cfg.block_size;
    
    // Redirect output to file
    std::ofstream ofs;
    std::streambuf* original_cout = std::cout.rdbuf();
    std::string output_file = folder + "/experiment_log.txt";
    ofs.open(output_file);
    if (ofs.is_open()) std::cout.rdbuf(ofs.rdbuf());
    
    // Logging header
    std::cout << "============================\n";
    std::cout << "Running " << pattern_str << " experiments\n";
    std::cout << "Matrix size: " << rows << " x " << cols << "\n";
    std::cout << "Block size: " << cfg.block_size << "\n";
    std::cout << "Causal: " << (cfg.apply_causal ? "yes" : "no") << "\n";
    if (cfg.pattern == MaskPattern::BigBird || cfg.pattern == MaskPattern::LongFormer)
        std::cout << "Sparsity: " << cfg.sparsity << "\n";
    if (cfg.pattern == MaskPattern::LongNet)
        std::cout << "Segment start: " << cfg.segment_start << "\n";
    if (cfg.matrix_file) std::cout << "Matrix file: " << cfg.matrix_file.value() << "\n";
    std::cout << "============================\n";
    
    // CSR arrays
    size_t* ia = nullptr;
    size_t* ja = nullptr;

    if (cfg.matrix_file) {
        // Load matrix from file
        mmdata* mm = (mmdata*) calloc(1, sizeof(mmdata));
        char* mm_filename = new char[cfg.matrix_file.value().size() + 1];
        std::strcpy(mm_filename, cfg.matrix_file.value().c_str());

        if (initialize_mm(mm_filename, mm) != 0) {
            std::cerr << "Error initializing matrix from file: " << mm_filename << "\n";
            free(mm_filename);
            return;
        }
        free(mm_filename);

        mm2csr(mm, ia, ja, rows, cols);
        freemm(mm);
    } else {
        // Generate mask based on pattern
        switch(cfg.pattern) {
            case MaskPattern::BigBird: {
                int blocks = static_cast<int>((cfg.sparsity * rows) / 8);
                int window_size = 2 * blocks; // Example default
                int num_global_tokens = 3 * blocks; // or user-defined
                int random_per_row = 3 * blocks;
                generate_bigbird_mask(rows, cols, window_size, num_global_tokens, random_per_row, cfg.apply_causal, ia, ja);
                break;
            }
            case MaskPattern::LongNet: {
                int window_size = 0; // Example local window
                std::vector<size_t> segment_sizes;
                std::vector<size_t> dilation_sizes;
                for (int seg = cfg.segment_start; seg < static_cast<int>(std::log2(rows)); ++seg) {
                    segment_sizes.push_back(static_cast<size_t>(std::pow(2, seg)));
                    dilation_sizes.push_back(static_cast<size_t>(std::pow(2, seg - cfg.segment_start)));
                }
                generate_longnet_mask(rows, cols, window_size, segment_sizes, dilation_sizes, cfg.apply_causal, ia, ja);
                break;
            }
            case MaskPattern::LongFormer: {
                int blocks = static_cast<int>(cfg.sparsity * rows);
                int window_size = 2 * blocks; // Local window
                int dilation = static_cast<int>(std::pow(2, rand() % (static_cast<int>(std::log2(rows)))));    // Example dilation
                int num_internal_global_tokens = 6 * blocks; // user-defined
                generate_longformer_mask(rows, cols, window_size, dilation, num_internal_global_tokens, cfg.apply_causal, ia, ja);
                break;
            }
        }
    }

    // Run BSR experiments
    bool RCM = false, HPRownet = false, HPSB = false, HPNBM = false, TwoConstraint = false;
    bool RCM_HPNBM = false, HPSB_HPNBM = false, HPRownet_HPNBM = true;

    run_experiments(RCM, HPRownet, HPSB, HPNBM, TwoConstraint,
                    RCM_HPNBM, HPSB_HPNBM, HPRownet_HPNBM,
                    ia, ja, rows, cfg.block_size);

    free(ia);
    free(ja);

    // Redirect output back
    std::cout.rdbuf(original_cout);
    if (ofs.is_open()) ofs.close();
}


//
// Help menu
//
void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [OPTIONS]\n\n"
              << "Required:\n"
              << "  -n <int>         Number of blocks\n"
              << "  -b <int>         Block size\n"
              << "  -p <pattern>     Mask pattern: bigbird | longnet | longformer\n\n"
              << "Optional:\n"
              << "  -s <double>      Sparsity (required for bigbird, longformer)\n"
              << "  -g <int>         Segment start exponent (required for longnet)\n"
              << "  --causal         Apply causal masking (default: off)\n"
              << "  -o <path>        Save output path (default: none)\n"
              << "  -h, --help       Show help\n\n"
              << "Examples:\n"
              << "  " << prog << " -n 4096 -b 128 -p bigbird -s 0.1 --causal\n"
              << "  " << prog << " -n 4096 -b 128 -p longnet -g 2 -o out/\n";
}

//
// Argument parser
//
bool parse_args(int argc, char* argv[], Config &cfg) {
    if (argc == 1) {
        print_usage(argv[0]);
        return false;
    }

    bool pattern_set = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return false;
        }
        else if (arg == "-n" && i + 1 < argc) {
            cfg.no_blocks = std::stoi(argv[++i]);
        }
        else if (arg == "-b" && i + 1 < argc) {
            cfg.block_size = std::stoi(argv[++i]);
        }
        else if (arg == "-p" && i + 1 < argc) {
            std::string p = argv[++i];
            if      (p == "bigbird")    { cfg.pattern = MaskPattern::BigBird;    pattern_set = true; }
            else if (p == "longnet")    { cfg.pattern = MaskPattern::LongNet;    pattern_set = true; }
            else if (p == "longformer") { cfg.pattern = MaskPattern::LongFormer; pattern_set = true; }
            else {
                std::cerr << "Unknown pattern: " << p << "\n";
                return false;
            }
        }
        else if (arg == "-s" && i + 1 < argc) {
            cfg.sparsity = std::stod(argv[++i]);
        }
        else if (arg == "-g" && i + 1 < argc) {
            cfg.segment_start = std::stoi(argv[++i]);
        }
        else if (arg == "--causal") {
            cfg.apply_causal = true;
        }
        else if (arg == "-o" && i + 1 < argc) {
            cfg.save_path = argv[++i];
        }
        else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            return false;
        }
    }

    //
    // Validation
    //
    if (cfg.no_blocks <= 0) {
        std::cerr << "Error: -n <no_blocks> is required and must be > 0\n";
        return false;
    }

    if (cfg.block_size <= 0) {
        std::cerr << "Error: -b <block_size> is required and must be > 0\n";
        return false;
    }

    if (!pattern_set) {
        std::cerr << "Error: -p <pattern> is required\n";
        return false;
    }

    // Pattern-specific argument requirements
    if (cfg.pattern == MaskPattern::BigBird ||
        cfg.pattern == MaskPattern::LongFormer)
    {
        if (cfg.sparsity <= 0) {
            std::cerr << "Error: -s <sparsity> is required for "
                      << pattern_to_string(cfg.pattern) << "\n";
            return false;
        }
    }

    if (cfg.pattern == MaskPattern::LongNet) {
        if (cfg.segment_start < 0) {
            std::cerr << "Error: -g <segment_start> is required for longnet\n";
            return false;
        }
    }

    return true;
}

int main(int argc, char* argv[]) {
    Config cfg;

    if (!parse_args(argc, argv, cfg)) {
        return 1;
    }
    // Run experiments
    run_experiments_mask(cfg);
    return 0; // Indicate successful execution
}