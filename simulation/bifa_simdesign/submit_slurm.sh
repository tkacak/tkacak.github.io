#!/bin/bash
#SBATCH --job-name=bifa_sim
#SBATCH --array=1-608
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/bifa_%A_%a.out
#SBATCH --error=logs/bifa_%A_%a.err

module load r

cd "$SLURM_SUBMIT_DIR"
export MAX_TIME="21:30:00"
export MAX_RAM="7GB"

Rscript run_hpc.R
