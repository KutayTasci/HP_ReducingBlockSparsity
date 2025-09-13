# HP_ReducingBlockSparsity

This repository provides implementations of **sparse matrix reordering techniques** based on the **Reverse Cuthill–McKee (RCM) algorithm** and various **hypergraph partitioning approaches**.  
The goal is to minimize the number of nonzero blocks in block-sparse matrix representations (e.g., BlockCSR), which can improve performance in scientific computing and high-performance linear algebra applications.

---

## ✨ Features

- **RCM (Reverse Cuthill–McKee)**  
  A classic bandwidth reduction method for sparse matrices.

- **HPSB (Hypergraph Partitioning for Sparse Blocks)**  
  Uses hypergraph partitioning to cluster rows/columns and reduce block fill-in.

- **HPNBM (Hypergraph Partitioning for Nonzero Blocks Minimization)**  
  Targets minimizing the number of nonzero blocks directly.

- **RCM + HPNBM Hybrid**  
  Applies RCM first, followed by block minimization.

- **HPSB + HPNBM Hybrid**  
  Combines hypergraph-based sparse block clustering with block minimization.

---

## 🛠️ Code Structure

The main entry point demonstrates how to use the different reordering strategies:

```cpp
BlockCSR* bcsr = nullptr;

// Reverse Cuthill–McKee
reorder_RCM(ia, ja, rows, block_size, bcsr, true);

// Hypergraph Partitioning for Sparse Blocks
bcsr = nullptr;
reorder_HPSB(ia, ja, rows, block_size, bcsr, true);

// Hypergraph Partitioning for Nonzero Blocks Minimization
bcsr = nullptr;
reorder_HPNBM(ia, ja, rows, block_size, bcsr, true);

// Hybrid: RCM + HPNBM
bcsr = nullptr;
reorder_RCM_HPNBM(ia, ja, rows, block_size, bcsr, true);

// Hybrid: HPSB + HPNBM
bcsr = nullptr;
reorder_HPSB_HPNBM(ia, ja, rows, block_size, bcsr, true);
