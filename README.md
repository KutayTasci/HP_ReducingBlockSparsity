# HP_ReducingBlockSparsity

This repository provides implementations of **sparse matrix reordering and block sparsity reduction techniques**, with a focus on **mask generation** and **hypergraph-based reordering**. The goal is to minimize the number of nonzero blocks in block-sparse matrices, improving performance in **high-performance computing** and **scientific computing applications**.

---

## 📦 Features

- **RCM (Reverse Cuthill–McKee)** — Classic bandwidth reduction for sparse matrices.  
- **HPRownet / HPSB / HPNBM** — Hypergraph partitioning methods for block sparsity reduction.  
- **Hybrid Approaches** — Combine RCM and hypergraph minimization strategies (e.g., `RCM + HPNBM`, `HPSB + HPNBM`, `HPRownet + HPNBM`).  
- **Mask Patterns for Sparse Matrices** — Supports structured sparsity patterns:
  - `bigbird` — Global + local windowed sparsity  
  - `longnet` — Dilated sliding window with segment start  
  - `longformer` — Dilated sliding window + internal global connections  

---

## 🛠️ Dependencies

### General
- **C++17** or newer  
- **CMake**  
- **Boost**  
- **Intel TBB (Threading Building Blocks)**  
- **hwloc** (Hardware locality library)  

### HPC-specific (Case Western HPC)
- `module load Mt_KaHyPar/hp`  
- `module load PaToH/hp`  
- `module load CMake/3.26.3-GCCcore-12.3.0`  
- `module load Boost/1.82.0-GCC-12.3.0`  
- `module load tbb/2021.11.0-GCCcore-12.3.0`  
- `module load hwloc/2.9.1-GCCcore-12.3.0`  

> ⚠️ The general dependencies (Boost, TBB, hwloc) are required on any system.

---

## 🏗️ Compilation

```bash
mkdir build
cd build
cmake ..
make
```

---

## 🚀 Usage
```
Usage: ./TestHP [OPTIONS]

Required:
  -n <int>         Number of blocks
  -b <int>         Block size
  -p <pattern>     Mask pattern: bigbird | longnet | longformer

Optional:
  -s <double>      Sparsity (required for bigbird, longformer)
  -g <int>         Segment start exponent (required for longnet)
  --causal         Apply causal masking (default: off)
  -o <path>        Save output path (default: none)
  -h, --help       Show help

Examples:
  ./TestHP -n 4096 -b 128 -p bigbird -s 0.1 --causal
  ./TestHP -n 4096 -b 128 -p longnet -g 2 -o out/

### Examples

# BigBird pattern with causal masking
./TestHP -n 4096 -b 128 -p bigbird -s 0.1 --causal

# LongNet pattern with segment start
./TestHP -n 4096 -b 128 -p longnet -g 2 -o out/

# LongFormer pattern with sparsity
./TestHP -n 4096 -b 128 -p longformer -s 0.05
```
---

## 📖 Notes

- Ensure that hypergraph partitioning libraries (mt-KaHyPar, PaToH) are installed if you plan to use HPNBM-based methods.
- Mask patterns (bigbird, longnet, longformer) are intended for structured sparse operations, commonly used in attention mechanisms in machine learning.
