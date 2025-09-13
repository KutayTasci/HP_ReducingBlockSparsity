#include "HPNBM/hp.h"
#include <mtkahypar.h>

#include <memory>
#include <vector>
#include <iostream>
#include <thread>
#include <cassert>



/**
 * @brief Partitions the columns of a sparse matrix using row-net hypergraph partitioning.
 *
 * This function constructs a row-net hypergraph from the given CSR (Compressed Sparse Row) representation
 * of a sparse matrix and partitions its columns into the specified number of parts using the mt-kahypar library.
 * The partitioning aims to minimize the number of cut hyperedges (rows connecting columns in different partitions).
 * The output col_part array contains the column indices, ordered such that columns with incident cut hyperedges
 * are placed at the end.
 *
 * @param[in] ia         Row pointer array of size (rows + 1), CSR format.
 * @param[in] ja         Column indices array, CSR format.
 * @param[in] rows       Number of rows in the matrix.
 * @param[in] cols       Number of columns in the matrix.
 * @param[in] num_parts  Desired number of partitions.
 * @param[out] col_part  Pointer to array of size 'cols' that will contain the reordered column indices after partitioning.
 *
 * @note
 * - Uses mt-kahypar for hypergraph partitioning.
 * - Initializes all node and net weights to 1.
 * - Checks for valid input arrays and bounds.
 * - Allocates and manages memory for intermediate arrays.
 * - The output col_part array must be deleted by the caller to avoid memory leaks.
 * - The function prints status messages and error information to stdout/stderr.
 *
 * @warning
 * - The function assumes that ia and ja are valid CSR arrays.
 * - If input arrays are invalid or out-of-bounds indices are detected, the function prints an error and returns early.
 * - The function may terminate the program if mt-kahypar fails to construct or partition the hypergraph.
 */
void HPSB_RowNet(
    size_t* ia,
    size_t* ja,
    size_t rows,
    size_t cols,
    size_t num_parts,
    size_t*& col_part
) {

    if (cols == 0 || ia[rows] == 0) {
        std::cerr << "Error: Zero-sized allocation detected." << std::endl;
        return;
    }

    mt_kahypar_error_t error{};
    col_part = new size_t[cols];

    // Initialize mt-kahypar
    mt_kahypar_initialize(std::thread::hardware_concurrency(), true);

    // Hypergraph construction
    const mt_kahypar_hypernode_id_t num_nodes = cols;
    const mt_kahypar_hyperedge_id_t num_hyperedges = rows;

    std::unique_ptr<size_t[]> net_indices = std::make_unique<size_t[]>(cols + 1);
    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> nets = std::make_unique<mt_kahypar_hyperedge_id_t[]>(ia[rows]);

    // Build net_indices
    std::fill(net_indices.get(), net_indices.get() + cols + 1, 0);
    for (size_t i = 0; i < ia[rows]; ++i) {
        if (ja[i] >= cols) {
            std::cerr << "Error: ja[" << i << "] = " << ja[i] << " is out of bounds (cols = " << cols << ")." << std::endl;
            return;
        }
        net_indices[ja[i] + 1]++;
    }
    for (size_t i = 1; i <= cols; ++i) {
        net_indices[i] += net_indices[i - 1];
    }

    // Build nets
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = ia[i]; j < ia[i + 1]; ++j) {
            if (ja[j] >= cols) {
                std::cerr << "Error: ja[" << j << "] = " << ja[j] << " is out of bounds (cols = " << cols << ")." << std::endl;
                return;
            }
            nets[net_indices[ja[j]]] = i;
            net_indices[ja[j]]++;
        }
    }
    for (size_t i = cols; i > 0; --i) {
        net_indices[i] = net_indices[i - 1];
    }
    net_indices[0] = 0;

    // Node and net weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> node_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(num_nodes);
    std::fill(node_weights.get(), node_weights.get() + num_nodes, 1);

    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> net_weights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(num_hyperedges);
    std::fill(net_weights.get(), net_weights.get() + num_hyperedges, 1);

    // Create context and hypergraph
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
        context, num_nodes, num_hyperedges,
        net_indices.get(), nets.get(),
        net_weights.get(), node_weights.get(), &error
    );
    if (!hypergraph.hypergraph) {
        std::cout << error.msg << std::endl;
        std::exit(1);
    }

    // Partitioning parameters
    mt_kahypar_set_partitioning_parameters(context, num_parts, 0.10, KM1);
    mt_kahypar_set_seed(42);
    mt_kahypar_status_t status = mt_kahypar_set_context_parameter(context, VERBOSE, "0", &error);
    assert(status == SUCCESS);

    // Partition hypergraph
    mt_kahypar_partitioned_hypergraph_t partitioned_hg = mt_kahypar_partition(hypergraph, context, &error);
    if (!partitioned_hg.partitioned_hg) {
        std::cout << error.msg << std::endl;
        std::exit(1);
    }

    // Extract partition
    std::unique_ptr<mt_kahypar_partition_id_t[]> partition = std::make_unique<mt_kahypar_partition_id_t[]>(num_nodes);
    mt_kahypar_get_partition(partitioned_hg, partition.get());

    // Compute part sizes
    std::unique_ptr<size_t[]> part_sizes = std::make_unique<size_t[]>(num_parts + 1);
    std::fill(part_sizes.get(), part_sizes.get() + num_parts + 1, 0);
    for (size_t i = 0; i < cols; ++i) {
        part_sizes[partition[i] + 1]++;
    }
    for (size_t i = 1; i <= num_parts; ++i) {
        part_sizes[i] += part_sizes[i - 1];
    }

    // Reorder columns by partition
    std::unique_ptr<size_t[]> temp_part = std::make_unique<size_t[]>(cols);
    for (size_t i = 0; i < cols; ++i) {
        int part = partition[i];
        temp_part[part_sizes[part]] = i;
        part_sizes[part]++;
    }

    // Place columns with cut hyperedges at the end
    size_t start = 0, end = cols - 1;
    for (size_t i = 0; i < cols; ++i) {
        mt_kahypar_hyperedge_id_t cut_count = mt_kahypar_num_incident_cut_hyperedges(partitioned_hg, temp_part[i]);
        if (cut_count > 0) {
            col_part[end--] = temp_part[i];
        } else {
            col_part[start++] = temp_part[i];
        }
    }

    // Cleanup
    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
}

/**
 * @brief Partitions the rows of a sparse matrix into blocks using hypergraph partitioning.
 *
 * This function constructs a custom hypergraph based on the block structure and partitions the rows
 * into the specified number of blocks using the mt-kahypar library. It supports fixed vertices for block anchors.
 * The output row_part array contains the reordered row indices after partitioning, and block_row_ptr
 * contains the starting indices of each block in row_part.
 *
 * @param[in] ia           Row pointer array of size (rows + 1), CSR format.
 * @param[in] ja           Column indices array, CSR format.
 * @param[in] rows         Number of rows in the matrix.
 * @param[in] cols         Number of columns in the matrix.
 * @param[in] num_parts    Desired number of partitions (blocks).
 * @param[in] block_size   Size of each block.
 * @param[out] row_part    Pointer to array of size 'rows' that will contain the reordered row indices.
 * @param[out] block_row_ptr Pointer to array of size (num_parts + 1) with block boundaries in row_part.
 *
 * @note
 * - Uses mt-kahypar for hypergraph partitioning.
 * - Initializes all node and net weights to 1, except block anchor nodes.
 * - Handles fixed vertices for block anchors.
 * - The output arrays must be deleted by the caller to avoid memory leaks.
 * - Prints status messages and error information to stdout/stderr.
 */
void HPNBM(
    size_t* ia,
    size_t* ja,
    size_t rows,
    size_t cols,
    size_t num_parts,
    size_t block_size,
    size_t*& row_part,
    size_t*& block_row_ptr
) {
    // Error checking for input arrays
    if (ia == nullptr || ja == nullptr || ia[rows] == 0) {
        std::cerr << "Error: Invalid input arrays or no non-zero entries." << std::endl;
        return;
    }

    mt_kahypar_error_t error{};
    row_part = new size_t[rows];

    // Initialize mt-kahypar with all available cores and NUMA policy
    mt_kahypar_initialize(std::thread::hardware_concurrency(), true);

    // Initial partition assignment for each row
    std::vector<size_t> initial_parts(rows, 0);
    for (size_t i = 0; i < rows; i++) {
        initial_parts[i] = i / block_size;
    }

    // Prepare data structures for hypergraph construction
    std::vector<size_t> nodes_map; // Maps rows to hypergraph nodes
    std::vector<std::vector<size_t>> nodes_part_map(num_parts); // Maps nodes to partitions
    std::vector<int> part_local_sizes(num_parts, 1); // Local sizes for block anchor nodes

    size_t num_nodes_t = 0;      // Number of non-anchor nodes
    size_t num_connections = 0;  // Number of connections in the hypergraph

    // Build hypergraph node mapping and partition connections
    for (size_t i = 0; i < rows; i++) {
        bool flag = false;
        std::vector<size_t> temp(num_parts, 0);

        for (size_t j = ia[i]; j < ia[i + 1]; j++) {
            size_t col = ja[j];
            if (initial_parts[i] != initial_parts[col]) {
                flag = true;
                temp[initial_parts[col]] = 1;
            } else {
                temp[initial_parts[i]] = 1;
            }
        }
        if (flag) {
            nodes_map.push_back(i);
            for (size_t k = 0; k < num_parts; k++) {
                if (temp[k] == 1) {
                    num_connections++;
                    nodes_part_map[k].push_back(num_nodes_t);
                }
            }
            num_nodes_t++;
        } else {
            part_local_sizes[initial_parts[i]]++;
        }
    }

    // Hypergraph parameters
    const mt_kahypar_hypernode_id_t num_nodes = num_nodes_t + num_parts; // Non-anchor nodes + anchor nodes
    const mt_kahypar_hyperedge_id_t num_hyperedges = num_parts;

    // Allocate net_indices and nets arrays for hypergraph
    std::unique_ptr<size_t[]> net_indices = std::make_unique<size_t[]>(num_parts + 1);
    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> nets = std::make_unique<mt_kahypar_hyperedge_id_t[]>(num_connections + num_parts);

    net_indices[0] = 0;
    for (size_t i = 0; i < num_parts; i++) {
        net_indices[i + 1] = net_indices[i] + nodes_part_map[i].size() + 1;
        for (size_t j = 0; j < nodes_part_map[i].size(); j++) {
            nets[net_indices[i] + j] = nodes_part_map[i][j];
        }
        // Add anchor node for each part
        nets[net_indices[i] + nodes_part_map[i].size()] = num_nodes_t + i;
    }

    // Debug output
    std::cout << "Hypergraph constructed.\n";
    std::cout << "num_nodes: " << num_nodes << ", num_hyperedges: " << num_hyperedges << std::endl;

    if (num_nodes == 0 || num_hyperedges == 0) {
        std::cerr << "Error: num_nodes or num_hyperedges is zero, cannot allocate weights arrays." << std::endl;
        return;
    }

    // Node weights: 1 for normal nodes, part_local_sizes for anchor nodes
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> node_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(num_nodes);
    for (size_t i = 0; i < num_nodes_t; i++) {
        node_weights[i] = 1;
    }
    for (size_t i = 0; i < num_parts; i++) {
        node_weights[num_nodes_t + i] = part_local_sizes[i];
    }

    // Net weights: all 1
    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> net_weights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(num_hyperedges);
    std::fill(net_weights.get(), net_weights.get() + num_hyperedges, 1);

    // Create mt-kahypar context and hypergraph
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(
        context, num_nodes, num_hyperedges,
        net_indices.get(), nets.get(),
        net_weights.get(), node_weights.get(), &error
    );
    if (!hypergraph.hypergraph) {
        std::cout << error.msg << std::endl;
        std::exit(1);
    }

    // Set fixed vertices for anchor nodes
    std::unique_ptr<mt_kahypar_partition_id_t[]> fixed_vertices = std::make_unique<mt_kahypar_partition_id_t[]>(num_nodes);
    for (size_t i = 0; i < num_nodes_t; i++) {
        fixed_vertices[i] = -1;
    }
    for (size_t i = 0; i < num_parts; i++) {
        fixed_vertices[num_nodes_t + i] = i;
    }
    mt_kahypar_status_t status = mt_kahypar_add_fixed_vertices(hypergraph, fixed_vertices.get(), num_parts, &error);
    if (status != SUCCESS) {
        std::cout << error.msg << std::endl;
        std::exit(1);
    }
    std::cout << "Fixed vertices set.\n";

    // Set partitioning parameters
    mt_kahypar_set_partitioning_parameters(context, num_parts, 0.1, KM1);
    mt_kahypar_set_seed(42);
    status = mt_kahypar_set_context_parameter(context, VERBOSE, "0", &error);
    assert(status == SUCCESS);

    // Partition the hypergraph
    mt_kahypar_partitioned_hypergraph_t partitioned_hg = mt_kahypar_partition(hypergraph, context, &error);
    if (!partitioned_hg.partitioned_hg) {
        std::cout << error.msg << std::endl;
        std::exit(1);
    }

    // Extract partition assignment
    std::unique_ptr<mt_kahypar_partition_id_t[]> partition = std::make_unique<mt_kahypar_partition_id_t[]>(num_nodes);
    mt_kahypar_get_partition(partitioned_hg, partition.get());

    // Compute metrics (imbalance, km1)
    const double imbalance = mt_kahypar_imbalance(partitioned_hg, context);
    const int km1 = mt_kahypar_km1(partitioned_hg);

    // Allocate part_sizes array for block boundaries
    size_t* part_sizes = new size_t[num_parts + 1];
    std::fill(part_sizes, part_sizes + num_parts + 1, 0);

    // Map partitioning back to original rows
    size_t* tmp_parts = new size_t[rows];
    size_t tmp = 0;
    for (size_t i = 0; i < rows; i++) {
        if (tmp < nodes_map.size() && i == nodes_map[tmp]) {
            tmp_parts[i] = partition[tmp];
            tmp++;
        } else {
            tmp_parts[i] = initial_parts[i];
        }
        part_sizes[tmp_parts[i] + 1]++;
    }
    if (tmp != num_nodes_t) {
        std::cerr << "Error: Mismatch in number of mapped nodes. Expected " << num_nodes_t << ", but got " << tmp << "." << std::endl;
        exit(EXIT_FAILURE);
    }
    for (size_t i = 1; i <= num_parts; i++) {
        part_sizes[i] += part_sizes[i - 1];
    }

    // Set block_row_ptr boundaries
    for (size_t i = 0; i <= num_parts; i++) {
        block_row_ptr[i] = part_sizes[i];
    }

    // Reorder rows by partition
    for (size_t i = 0; i < rows; i++) {
        int part = tmp_parts[i];
        row_part[part_sizes[part]] = i;
        part_sizes[part]++;
    }

    // Cleanup
    delete[] part_sizes;
    delete[] tmp_parts;
    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
}
