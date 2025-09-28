#!/bin/bash
#SBATCH -p batch
#SBATCH -A sxk1942
#SBATCH -c 32
#SBATCH --mem=156G
#SBATCH -t 24:00:00
#SBATCH -J my_job
#SBATCH -o my_job_%j.out
#SBATCH -e my_job_%j.err

echo "Starting job on `date`"
echo "Running on node `hostname`"
echo "Current directory is `pwd`"

ulimit -c unlimited

#/home/kxt437/BlockMin/HP_ReducingBlockSparsity/build/TestHP 4096 0.0 0
#/home/kxt437/BlockMin/HP_ReducingBlockSparsity/build/TestHP 4096 0.0 0
/home/kxt437/BlockMin/HP_ReducingBlockSparsity/build/TestHP 4096 0.0 10
/home/kxt437/BlockMin/HP_ReducingBlockSparsity/build/TestHP 4096 0.0 11

