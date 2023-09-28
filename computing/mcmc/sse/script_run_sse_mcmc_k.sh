#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3850
#SBATCH --time=20:00:00
#SBATCH --partition=stats
#SBATCH --array=1
srun Rscript SCRIPT_RUN_SSE_MCMC_K_V2.R $SLURM_ARRAY_TASK_ID