#!/bin/bash
#SBATCH -p batch
#SBATCH -A sxk1942
#SBATCH -c 32
#SBATCH --mem=240G
#SBATCH -t 24:00:00
#SBATCH -J my_job
#SBATCH -o my_job_%j.out
#SBATCH -e my_job_%j.err

echo "Starting job on $(date)"
echo "Running on node $(hostname)"
echo "Current directory is $(pwd)"

ulimit -c unlimited

# Define paths
EXE_SRC="/home/kxt437/BlockMin/HP_ReducingBlockSparsity/build/TestHP_Long"
EXE_COPY="./TestHP_${SLURM_JOB_ID}"

# Copy the executable with unique name
cp "$EXE_SRC" "$EXE_COPY"
chmod +x "$EXE_COPY"

# Run experiments with the copied executable
"$EXE_COPY" 2048 0.0001 8
"$EXE_COPY" 2048 0.0005 9
"$EXE_COPY" 2048 0.001 10
"$EXE_COPY" 2048 0.005 11

# Remove the copied executable after runs
rm -f "$EXE_COPY"

echo "Job finished on $(date)"
