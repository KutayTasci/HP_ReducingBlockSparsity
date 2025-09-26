#include "HPNBM/util.h"


#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <iostream>
//#include <omp.h>


struct localBlock {
    size_t id;
    size_t row;
    size_t col;
    size_t nnz;
    size_t local_rows;
    std::vector<size_t> ia_coo;
    std::vector<size_t> ja_coo;
};


void testBlockSparsity(size_t* block_row_ptr, size_t* block_col_ind, size_t num_blocks, size_t* ia, size_t* ja, size_t rows, size_t cols) {

    size_t total_blocks = num_blocks * num_blocks;
    size_t non_empty_blocks = 0;
    size_t i, j;

    size_t* blocks = new size_t[total_blocks]();

    for (i = 0; i < total_blocks; ++i) {
        blocks[i] = 0;
    }

    size_t row_block = 0; size_t col_block = 0;
    for (i = 0; i < rows; i++) {
        if (i >= block_row_ptr[row_block + 1]) {
            row_block++;
        }
        
        for (j = ia[i]; j < ia[i + 1]; ++j) {
            size_t col = ja[j];
            col_block = 0;
            while (col >= block_col_ind[col_block + 1]) {
                col_block++;
            }
            size_t block_index = row_block * num_blocks + col_block;
            if (blocks[block_index] == 0) {
                blocks[block_index] = 1;
                non_empty_blocks++;
            }
        }
    }

    double sparsity = (double)non_empty_blocks / total_blocks;
    std::cout << "Total blocks: " << total_blocks << std::endl;
    std::cout << "Non-empty blocks: " << non_empty_blocks << std::endl;
    std::cout << "Block sparsity: " << sparsity * 100 << "%" << std::endl;
}

// Reverse Cuthill-McKee algorithm implementation with symmetric treatment
void ReverseCuthillMcKee(size_t* ia, size_t* ja, size_t n, size_t*& perm) {
    // Build undirected adjacency list
    
    std::vector<std::vector<size_t>> adj(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t idx = ia[i]; idx < ia[i + 1]; ++idx) {
            size_t j = ja[idx];
            adj[i].push_back(j);
            if (j < n) {
                adj[j].push_back(i); // Add reverse edge for undirected graph
            }
        }
    }

    // Remove duplicate neighbors
    for (size_t i = 0; i < n; ++i) {
        std::sort(adj[i].begin(), adj[i].end());
        adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
    }
    
    std::vector<bool> visited(n, false);
    std::vector<size_t> degrees(n, 0);
    for (size_t i = 0; i < n; ++i) {
        degrees[i] = adj[i].size();
    }

    std::vector<size_t> rcm_order;
    rcm_order.reserve(n);

    
    bool stop = false;
    // Handle disconnected components
    while (stop == false) {
        // find the unvisited node with the smallest degree
        size_t start = n; // Invalid index
        size_t min_degree = n + 1; // Larger than any possible degree
        for (size_t i = 0; i < n; ++i) {
            if (!visited[i] && degrees[i] < min_degree) {
                min_degree = degrees[i];
                start = i;
            }
        }
        if (start == n) {
            stop = true;
            continue; // All nodes visited
        }
        if (visited[start]) continue;

        std::queue<size_t> q;
        q.push(start);
        visited[start] = true;

        while (!q.empty()) {
            size_t v = q.front();
            q.pop();
            rcm_order.push_back(v);

            // Gather unvisited neighbors
            std::vector<size_t> neighbors;
            for (size_t u : adj[v]) {
                if (!visited[u]) {
                    neighbors.push_back(u);
                    visited[u] = true;
                }
            }
            // Sort neighbors by degree
            std::sort(neighbors.begin(), neighbors.end(),
                      [&](size_t a, size_t b) { return degrees[a] < degrees[b]; });

            for (size_t u : neighbors) {
                q.push(u);
            }
        }
    }
    
    // Reverse the order for RCM
    perm = new size_t[n];
    for (size_t i = 0; i < n; ++i) {
        perm[i] = rcm_order[n - i - 1];
    }

    // Permutation shows new index for each old index
    // i.e., old index perm[i] goes to new index i

    // No need to manually delete std::vector memory; vectors are automatically managed.
    // clean up adjacency list
}

void reorderCSR(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* row_perm, size_t* col_perm, size_t*& ia_new, size_t*& ja_new) {
    ia_new = new size_t[rows + 1];
    ia_new[0] = 0;
    size_t nnz = ia[rows];
    ja_new = new size_t[nnz];

    // Create a mapping from old column indices to new column indices
    std::vector<size_t> col_map(cols);
    for (size_t j = 0; j < cols; ++j) {
        col_map[col_perm[j]] = j;
    }

    for (size_t i = 0; i < rows; ++i) {
        size_t old_row = row_perm[i];
        size_t start = ia[old_row];
        size_t end = ia[old_row + 1];
        ia_new[i + 1] = ia_new[i] + (end - start);

        for (size_t idx = start; idx < end; ++idx) {
            size_t old_col = ja[idx];
            size_t new_col = col_map[old_col];
            ja_new[ia_new[i] + (idx - start)] = new_col;
        }
    }

}

void reorderCSR_Rows(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* row_perm, size_t*& ia_new, size_t*& ja_new) {
    ia_new = new size_t[rows + 1];
    ia_new[0] = 0;
    size_t nnz = ia[rows];
    ja_new = new size_t[nnz];

    for (size_t i = 0; i < rows; ++i) {
        size_t old_row = row_perm[i];
        if (old_row >= rows) {
            std::cerr << "Error: row_perm[" << i << "] = " << old_row << " is out of bounds (rows = " << rows << ")." << std::endl;
            exit(EXIT_FAILURE);
        }
        size_t start = ia[old_row];
        size_t end = ia[old_row + 1];
        ia_new[i + 1] = ia_new[i] + (end - start);

        for (size_t idx = start; idx < end; ++idx) {
            ja_new[ia_new[i] + (idx - start)] = ja[idx];
        }
    }
}

void reorderCSR_Cols(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* col_perm, size_t*& ia_new, size_t*& ja_new) {
    // just copy ia as is
    ia_new = new size_t[rows + 1];
    for (size_t i = 0; i <= rows; ++i) {
        ia_new[i] = ia[i];
    }
    
    size_t nnz = ia[rows];
    ja_new = new size_t[nnz];
    
    // Create a mapping from old column indices to new column indices
    std::vector<size_t> col_map(cols);
    for (size_t j = 0; j < cols; ++j) {
        if (col_perm[j] >= cols) {
            std::cerr << "Error: col_perm[" << j << "] = " << col_perm[j] << " is out of bounds (cols = " << cols << ")." << std::endl;
            exit(EXIT_FAILURE);
        }
        col_map[col_perm[j]] = j;
    }
    //printf("Column reordering done.\n");
    for (size_t i = 0; i < ia[rows]; ++i) {
        size_t old_col = ja[i];
        size_t new_col = col_map[old_col];
        ja_new[i] = new_col;
    }
    
}

void create_BlockCSR(size_t* block_row_ind, size_t* block_col_ind, size_t num_blocks, size_t* ia, size_t* ja, size_t rows, size_t cols, size_t* row_perm, size_t* column_perm, BlockCSR*& bcsr) {
    size_t i, j;

    size_t* row_part = new size_t[rows];
    size_t* col_part = new size_t[cols];
    size_t block_row = 0;
    for (i = 0; i < rows; ++i) {
        size_t r = row_perm[i];
        
        while (i >= block_row_ind[block_row + 1]) {
            block_row++;
        }
        row_part[r] = block_row; // map original row index to block
    }

    size_t block_col = 0;
    for (j = 0; j < cols; ++j) {
        size_t c = column_perm[j];
        
        while (j >= block_col_ind[block_col + 1]) {
            block_col++;
        }
        col_part[c] = block_col; // map original col index to block
    }

    std::unordered_map<int, localBlock> blocks;

    size_t num_blocks_found = 0;
    for (i = 0; i < rows; i++)
    {
        int block_row = row_part[i];
        for (j = ia[i]; j < ia[i + 1]; j++)
        {
            size_t c = ja[j];
            int block_col = col_part[c];
            int block_id = (block_row * num_blocks) + block_col;

            if (blocks.find(block_id) == blocks.end()) {
                localBlock lb;
                lb.id = block_id;
                lb.row = block_row;
                lb.col = block_col;
                lb.nnz = 0;
                lb.local_rows = 0;
                lb.ia_coo.clear();
                lb.ja_coo.clear();
                blocks[block_id] = lb;
                num_blocks_found++;
            }
            blocks[block_id].ia_coo.push_back(i);
            blocks[block_id].ja_coo.push_back(c);
            blocks[block_id].nnz++;
        }
    }

    std::vector<std::pair<int, localBlock>> sorted_blocks(blocks.begin(), blocks.end());
    std::sort(sorted_blocks.begin(), sorted_blocks.end(),
              [](const std::pair<int, localBlock>& a, const std::pair<int, localBlock>& b) {
                  return a.first < b.first;
              });
    
    size_t ctr_rows = 0;
    size_t ctr_nnz = 0;
    //fill local_rows for each block
    for (auto& pair : sorted_blocks) {
        localBlock& lb = pair.second;
        //get the number of unique rows in ia_coo without modifying the original vector
        std::vector<size_t> unique_rows = lb.ia_coo;
        std::sort(unique_rows.begin(), unique_rows.end());
        auto last = std::unique(unique_rows.begin(), unique_rows.end());
        unique_rows.erase(last, unique_rows.end());
        lb.local_rows = unique_rows.size();
        ctr_rows += lb.local_rows;
        ctr_nnz += lb.nnz;
    }
    bcsr = new BlockCSR;
    bcsr->num_block_rows = num_blocks;
    bcsr->num_block_cols = num_blocks;
    bcsr->num_blocks = num_blocks_found;
    bcsr->block_indices = new size_t[num_blocks_found + 1];
    bcsr->block_indices[0] = 0;
    bcsr->block_id_map = new size_t[num_blocks_found];
    bcsr->ia = new size_t[ctr_rows + 1];
    bcsr->ia[0] = 0;
    bcsr->ja = new size_t[ctr_nnz];
    bcsr->row_map = new size_t[ctr_rows];
    bcsr->block_row_ind = new size_t[num_blocks + 1];
    bcsr->block_col_ind = new size_t[num_blocks + 1];
    bcsr->row_perm = new size_t[rows];
    bcsr->col_perm = new size_t[cols];
    for (i = 0; i < rows; i++) {
        bcsr->row_perm[i] = row_perm[i];
    }
    for (j = 0; j < cols; j++) {
        bcsr->col_perm[j] = column_perm[j];
    }

    for (i = 0; i <= num_blocks; i++) {
        bcsr->block_row_ind[i] = block_row_ind[i];
        bcsr->block_col_ind[i] = block_col_ind[i];
    }

    size_t ctr_block = 0; // block counter
    for (auto& pair : sorted_blocks) {
        localBlock& lb = pair.second;
        size_t block_idx = lb.id;
        bcsr->block_id_map[ctr_block] = block_idx;
        bcsr->block_indices[ctr_block + 1] = bcsr->block_indices[ctr_block] + lb.local_rows;

        size_t current_index = -1; 
        size_t current_row = -1;

        for (i = 0; i < lb.nnz; i++) {
            size_t r = lb.ia_coo[i];
            size_t c = lb.ja_coo[i];
            if (r != current_row) {
                current_row = r;
                current_index++;
                bcsr->row_map[bcsr->block_indices[ctr_block] + current_index] = r;
                bcsr->ia[bcsr->block_indices[ctr_block] + current_index + 1] = bcsr->ia[bcsr->block_indices[ctr_block] + current_index];
            }
            
            
            bcsr->ja[bcsr->ia[bcsr->block_indices[ctr_block] + current_index + 1]++] = c;
        }

        ctr_block++;
    }

    //print bcsr details
    /*
    std::cout << "Block Indices: ";
    for (i = 0; i <= num_blocks_found; i++) {
        printf("%zu ", bcsr->block_indices[i]);
    }
    std::cout << std::endl;
    std::cout << "IA:";
    for (i = 0; i <= ctr_rows; i++) {
        printf("%zu ", bcsr->ia[i]);
    }
    std::cout << " " << std::endl;
    std::cout << "Row mapping of all rows in all blocks: ";
    for (i = 0; i < ctr_rows; i++) {
        printf("%zu ", bcsr->row_map[i]);
    }
    std::cout << std::endl;
    std::cout  << "JA: ";
    for (i = 0; i < ctr_nnz; i++) {
        printf("%zu ", bcsr->ja[i]);
    }
    std::cout << std::endl;
    */

}

void create_BSR(size_t* ia, size_t* ja, size_t rows, size_t cols, size_t block_size, BSR*& bsr) {
    size_t i, j, k;

    size_t num_block_rows = (rows + block_size - 1) / block_size;
    size_t num_block_cols = (cols + block_size - 1) / block_size;

    size_t* block_offsets = new size_t[num_block_rows + 1];
    size_t num_blocks = 0;
    block_offsets[0] = 0;
    std::vector<size_t> block_col_indices_map;
    for (i = 0; i < num_block_rows; ++i) {
        size_t row_start = i * block_size;
        size_t row_end = std::min(row_start + block_size, rows);
        std::unordered_map<size_t, bool> block_cols_set;

        for (j = row_start; j < row_end; ++j) {
            for (size_t idx = ia[j]; idx < ia[j + 1]; ++idx) {
                size_t col = ja[idx];
                size_t block_col = col / block_size;
                block_cols_set[block_col] = true;
            }
        }
        num_blocks += block_cols_set.size();
        block_offsets[i + 1] = block_offsets[i] + block_cols_set.size();
        for (const auto& pair : block_cols_set) {
            block_col_indices_map.push_back(pair.first);
        }
        
    }
    size_t* block_col_indices = block_col_indices_map.data();

    int* row_major_values = new int[num_blocks * block_size * block_size];
    for (i = 0; i < num_blocks * block_size * block_size; ++i) {
        row_major_values[i] = 0;
    }

    // Fill in the BSR values with actual non-zero values from the original matrix
    
    for (i = 0; i < num_block_rows; ++i) {
        size_t row_start = i * block_size;
        size_t row_end = std::min(row_start + block_size, rows);
        std::unordered_map<size_t, std::vector<size_t>> block_col_to_index;
        for (j = row_start; j < row_end; ++j) {
            for (size_t idx = ia[j]; idx < ia[j + 1]; ++idx) {
                size_t col = ja[idx];
                size_t block_col = col / block_size;
                size_t local_row = j - row_start;
                size_t local_col = col % block_size;
                size_t index = local_row * block_size + local_col;
                block_col_to_index[block_col].push_back(index);
            }
        }

        size_t block_start = block_offsets[i];
        size_t block_end = block_offsets[i + 1];
        for (j = block_start; j < block_end; ++j) {
            size_t block_col = block_col_indices[j];
            for (const auto& index : block_col_to_index[block_col]) {
                row_major_values[j * block_size * block_size + index] = 1; // Assuming non-zero value is 1
            }
        }

    }
    bsr = new BSR;
    bsr->block_offsets = block_offsets;
    bsr->block_col_ind = block_col_indices;
    bsr->block_size = block_size;
    bsr->row_major_values = row_major_values;
    bsr->rows = rows;
    bsr->cols = cols;
    bsr->nnz = ia[rows];
    bsr->num_blocks = num_blocks;
}

void analyzeBlockCSR(BlockCSR* bcsr) {

    size_t total_nnz = bcsr->ia[bcsr->block_indices[bcsr->num_blocks]];
    double avg_block_nnz = (double)total_nnz / bcsr->num_blocks;
    double block_sparsity = (double)bcsr->num_blocks / (bcsr->num_block_rows * bcsr->num_block_cols);

    std::unordered_map<std::string, int> block_row_dimension_counts;

    for (size_t b = 0; b < bcsr->num_block_rows; ++b) {
        size_t start_row = bcsr->block_row_ind[b];
        size_t end_row = bcsr->block_row_ind[b + 1];

        size_t start_col = bcsr->block_col_ind[b];
        size_t end_col = bcsr->block_col_ind[b + 1];
        
        size_t block_rows = end_row - start_row;
        size_t block_cols = end_col - start_col;

        std::string key = std::to_string(block_rows) + "x" + std::to_string(block_cols);
        block_row_dimension_counts[key]++;
    }

    int* nnz_per_block = new int[bcsr->num_blocks];
    double avg_nnz_per_block = 0.0;
    double variance_nnz_per_block = 0.0;
    int max_nnz_per_block = 0;
    int min_nnz_per_block = total_nnz;
    for (size_t b = 0; b < bcsr->num_blocks; ++b) {
        size_t start = bcsr->block_indices[b];
        size_t end = bcsr->block_indices[b + 1];
        nnz_per_block[b] = bcsr->ia[end] - bcsr->ia[start];

        avg_nnz_per_block += nnz_per_block[b];
        if (nnz_per_block[b] > max_nnz_per_block) {
            max_nnz_per_block = nnz_per_block[b];
        }
        if (nnz_per_block[b] < min_nnz_per_block) {
            min_nnz_per_block = nnz_per_block[b];
        }
        variance_nnz_per_block += nnz_per_block[b] * nnz_per_block[b];
    }
    avg_nnz_per_block /= bcsr->num_blocks;
    variance_nnz_per_block = (variance_nnz_per_block / bcsr->num_blocks) - (avg_nnz_per_block * avg_nnz_per_block);

    //std::cout << "----------------------------------------" << std::endl;
    std::cout << "Block CSR Analysis Results:" << std::endl;
    //std::cout << "Number of block rows: " << bcsr->num_block_rows << std::endl;
    //std::cout << "Number of block columns: " << bcsr->num_block_cols << std::endl;
    std::cout << "Number of non-empty blocks: " << bcsr->num_blocks << std::endl;
    //std::cout << "Total number of non-zeros in all blocks: " << total_nnz << std::endl;
    std::cout << "Average number of non-zeros per block: " << avg_nnz_per_block << std::endl;
    //std::cout << "Variance of non-zeros per block: " << variance_nnz_per_block << std::endl;
    //std::cout << "Max number of non-zeros in a block: " << max_nnz_per_block << std::endl;
    //std::cout << "Min number of non-zeros in a block: " << min_nnz_per_block << std::endl;


    std::cout << "Block sparsity (non-empty blocks / total blocks): " << block_sparsity * 100 << "%" << std::endl;
    for (const auto& pair : block_row_dimension_counts) {
        std::cout << "Block dimension " << pair.first << ": " << pair.second << " blocks" << std::endl;
    }
}

void analyzeBSR(BSR* bsr) {
    size_t total_blocks = bsr->num_blocks;
    size_t block_size = bsr->block_size;
    size_t total_nnz = bsr->nnz;

    double avg_nnz_per_block = (double)total_nnz / total_blocks;
    double block_density = (double)total_nnz / (total_blocks * block_size * block_size);

    //std::cout << "----------------------------------------" << std::endl;
    std::cout << "BSR Analysis Results:" << std::endl;
    std::cout << "Number of non-empty blocks: " << total_blocks << std::endl;
    std::cout << "Average number of non-zeros per block: " << avg_nnz_per_block << std::endl;
    std::cout << "Block sparsity (non-empty blocks / total blocks): " << (double)total_blocks / ((bsr->rows + block_size - 1) / block_size * (bsr->cols + block_size - 1) / block_size) * 100 << "%" << std::endl;

    std::cout << "Block dimension: " << block_size << "x" << block_size << ":" << "-1 blocks" << std::endl;
}

void freeBlockCSR(BlockCSR& bcsr) {
    if (bcsr.block_indices) {
        delete[] bcsr.block_indices;
        bcsr.block_indices = nullptr;
    }
    if (bcsr.block_id_map) {
        delete[] bcsr.block_id_map;
        bcsr.block_id_map = nullptr;
    }
    if (bcsr.ia) {
        delete[] bcsr.ia;
        bcsr.ia = nullptr;
    }
    if (bcsr.row_map) {
        delete[] bcsr.row_map;
        bcsr.row_map = nullptr;
    }
    if (bcsr.block_row_ind) {
        delete[] bcsr.block_row_ind;
        bcsr.block_row_ind = nullptr;
    }
    if (bcsr.block_col_ind) {
        delete[] bcsr.block_col_ind;
        bcsr.block_col_ind = nullptr;
    }
    if (bcsr.row_perm) {
        delete[] bcsr.row_perm;
        bcsr.row_perm = nullptr;
    }
    if (bcsr.col_perm) {
        delete[] bcsr.col_perm;
        bcsr.col_perm = nullptr;
    }
    if (bcsr.ja) {
        delete[] bcsr.ja;
        bcsr.ja = nullptr;
    }
    if (bcsr.row_map) {
        delete[] bcsr.row_map;
        bcsr.row_map = nullptr;
    }
}

void freeBSR(BSR*& bsr) {
    if (bsr->block_offsets) {
        delete[] bsr->block_offsets;
        bsr->block_offsets = nullptr;
    }
    // Do NOT delete[] bsr->block_col_ind, as it points to memory managed by a std::vector
    bsr->block_col_ind = nullptr;
    if (bsr->row_major_values) {
        delete[] bsr->row_major_values;
        bsr->row_major_values = nullptr;
    }
    delete bsr;
}