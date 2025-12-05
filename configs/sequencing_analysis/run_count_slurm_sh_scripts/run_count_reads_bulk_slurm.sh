#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --qos=cpu_qos
#SBATCH --mem=32GB
#SBATCH --time=0-10:00:00
#SBATCH --output=./run_logs/run_count_reads_bulk_%A.out
#SBATCH --error=./run_logs/run_count_reads_bulk_%A.err

mkdir -p ./logs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate py312-tcr-toolbox

tcr_toolbox count-reads-bulk $1
