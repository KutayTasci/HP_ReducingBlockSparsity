#include "HPNBM/mask_generator.h"

#include <vector>
#include <random>
#include <set>
#include <iostream>
#include <cmath>

void add_self_connections(size_t rows, size_t cols, std::set<std::pair<size_t, size_t>>& mask_set) {
    for (size_t i = 0; i < std::min(rows, cols); ++i) {
        std::pair<size_t, size_t> rc;
        rc.first = i;
        rc.second = i;
        mask_set.insert(rc);
    }
}

void external_global_connections(size_t rows, size_t cols, int num_connections, std::set<std::pair<size_t, size_t>>& mask_set) {
    if (num_connections <= 0) return;

    std::default_random_engine generator;
    std::uniform_int_distribution<size_t> row_distribution(0, rows - 1);
    std::uniform_int_distribution<size_t> col_distribution(0, cols - 1);

    for (int i = 0; i < num_connections; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::pair<size_t, size_t> rc;
            rc.first = i;
            rc.second = j;
            mask_set.insert(rc);
        }

        for (int j = 0; j < rows; ++j) {
            std::pair<size_t, size_t> rc;
            rc.first = j;
            rc.second = i;
            mask_set.insert(rc);
        }
    }
}

void internal_global_connections(size_t rows, size_t cols, int num_connections, std::set<std::pair<size_t, size_t>>& mask_set) {
    if (num_connections <= 0) return;

    std::default_random_engine generator;
    generator.seed(std::random_device()());
    std::uniform_int_distribution<size_t> row_distribution(0, rows - 1);


    for (int i = 0; i < num_connections; ++i) {
        size_t node = row_distribution(generator);

        for (int j = 0; j < cols; ++j) {
            std::pair<size_t, size_t> rc;
            rc.first = node;
            rc.second = j;
            mask_set.insert(rc);
        }

        for (int j = 0; j < rows; ++j) {
            std::pair<size_t, size_t> rc;
            rc.first = j;
            rc.second = node;
            mask_set.insert(rc);
        }
    }
}

void sliding_window_connections(size_t rows, size_t cols, int window_size, int dilation, std::set<std::pair<size_t, size_t>>& mask_set) {
    if (window_size <= 0) return;

    for (size_t i = 0; i < rows; ++i) {
        for (int j = -window_size; j <= window_size; ++j) {
            int col = static_cast<int>(i) + j * dilation;
            if (col >= 0 && col < static_cast<int>(cols)) {
                std::pair<size_t, size_t> rc;
                rc.first = i;
                rc.second = col;
                mask_set.insert(rc);
            }
        }
    }
}

void random_connections(size_t rows, size_t cols, int connections_per_row, std::set<std::pair<size_t, size_t>>& mask_set) {
    if (connections_per_row <= 0) return;
    //set random seed to real random
    std::default_random_engine generator;
    generator.seed(std::random_device()());
    std::uniform_int_distribution<size_t> col_distribution(0, cols - 1);

    for (size_t i = 0; i < rows; ++i) {
        std::set<size_t> selected_cols;
        while (selected_cols.size() < static_cast<size_t>(connections_per_row)) {
            // check if generated column is already selected
            size_t col = col_distribution(generator);
            if (selected_cols.find(col) != selected_cols.end()) {
                continue;
            }
            selected_cols.insert(col);
        }
        for (const auto& col : selected_cols) {
            std::pair<size_t, size_t> rc;
            rc.first = i;
            rc.second = col;
            mask_set.insert(rc);
        }
    }
}

void remove_causal_connections(size_t rows, size_t cols, std::set<std::pair<size_t, size_t>>& mask_set) {
    //iterate through the set and remove connections where col > row
    for (auto it = mask_set.begin(); it != mask_set.end(); ) {
        if (it->second > it->first) {
            it = mask_set.erase(it);
        } else {
            ++it;
        }
    }
}

void TwoD_dilation(size_t rows, size_t cols, size_t* segment_sizes, size_t* dilation_sizes, std::set<std::pair<size_t, size_t>>& mask_set) {
    if (segment_sizes == nullptr || dilation_sizes == nullptr) return;
    int size_segment = std::log2(rows) - 2;
    int size_dilation = std::log2(rows)- 2;

    if (rows != cols) return;
    int n = rows;

    if (size_segment != size_dilation) return;

    // iterate trough segment sizes
    for (int seg = 0; seg < size_segment; ++seg) {
        int segment_size = segment_sizes[seg];
        int dilation_size = dilation_sizes[seg];
        for (int ptr = 0; ptr < n; ptr += segment_size) {
            for (int i = ptr; i < std::min(ptr + segment_size, n); i += dilation_size) {
                for (int j = ptr; j < std::min(ptr + segment_size, n); j += dilation_size) {
                    std::pair<size_t, size_t> rc;
                    rc.first = i;
                    rc.second = j;
                    mask_set.insert(rc);
                }
            }
        }
    }
}

void convert_mask_to_csr(size_t rows, size_t cols, const std::set<std::pair<size_t, size_t>>& mask_set, size_t*& ia, size_t*& ja) {
    ia = new size_t[rows + 1];
    ia[0] = 0;
    ja = new size_t[mask_set.size()];

    size_t current_row = 0;
    size_t nnz_counter = 0;
    for (const auto& rc : mask_set) {
        //std::cout << "(" << rc.first << ", " << rc.second << ")\n";
        while (current_row < rc.first) {
            ia[current_row + 1] = nnz_counter;
            current_row++;
        }
        ja[nnz_counter++] = rc.second;
    }
    // fill in the rest of ia
    while (current_row < rows) {
        ia[current_row + 1] = nnz_counter;
        current_row++;
    }
}



void generate_mask(
    size_t rows, 
    size_t cols, 
    int external_global, 
    int internal_global, 
    int sliding_window, 
    int dilation, 
    int random_per_row, 
    size_t* segment_sizes, 
    size_t* dilation_sizes,
    bool apply_causal, 
    size_t*& ia, 
    size_t*& ja
) {


    std::set<std::pair<size_t, size_t>> mask_set;
    std::default_random_engine generator;
    std::uniform_int_distribution<size_t> distribution(0, cols - 1);

    add_self_connections(rows, cols, mask_set);
    
    if (external_global > 0) {
        external_global_connections(rows, cols, external_global, mask_set);
    }

    if (internal_global > 0) {
        internal_global_connections(rows, cols, internal_global, mask_set);
    }

    if (sliding_window > 0) {
        sliding_window_connections(rows, cols, sliding_window, dilation, mask_set);
    }

    if (random_per_row > 0) {
        random_connections(rows, cols, random_per_row, mask_set);
    }
    
    if (segment_sizes != nullptr && dilation_sizes != nullptr) {
        TwoD_dilation(rows, cols, segment_sizes, dilation_sizes, mask_set);
    }
    printf("Self connections added.\n");
    if (apply_causal) {
        remove_causal_connections(rows, cols, mask_set);
    }
    

    convert_mask_to_csr(rows, cols, mask_set, ia, ja);
    std::cout << "Total non-zeros in the mask: " << mask_set.size() << std::endl;
    std::cout << "Sparsity of the mask: " << static_cast<double>(mask_set.size()) / (rows * cols) << std::endl;

}