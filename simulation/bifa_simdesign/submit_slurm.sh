#!/bin/bash
## Gonderim: mkdir -p logs && sbatch submit_slurm.sh
## (SLURM logs/ klasorunu kendisi olusturmaz; yoksa gorev sessizce duser.)
## --array ust siniri = ceiling(8748 / ROWS_PER_JOB). ROWS_PER_JOB=12 -> 729.
## Kumenizin MaxArraySize degeri daha kucukse ROWS_PER_JOB'u buyutun.
#SBATCH --job-name=bifa_sim
#SBATCH --array=1-729
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/bifa_%A_%a.out
#SBATCH --error=logs/bifa_%A_%a.err

module load r        # kumenize gore duzenleyin (orn. module load R/4.4.1)

export SIM_DIR="$SLURM_SUBMIT_DIR"
export OUT_DIR="${SCRATCH:-$SLURM_SUBMIT_DIR}/bifa_sim_results"
export ROWS_PER_JOB=12
export MAX_TIME="21:30:00"   # --time'in ~%90'i: sure dolmadan tamamlanan kosullar kaydedilir
export MAX_RAM="7GB"         # --mem'in ~%90'i

Rscript "$SIM_DIR/run_hpc.R"
