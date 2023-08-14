#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=3850
#SBATCH --time=4:00:00
#SBATCH --partition=stats
#SBATCH --array=1
srun Rscript SCRIPT_RUN_SSE_MCMC.R $SLURM_ARRAY_TASK_ID