#!/bin/bash 
#SBATCH --nodes=1                        # requests 6 compute servers
#SBATCH --ntasks-per-node=1              # runs 6 tasks on each server
#SBATCH --cpus-per-task=10              # uses 1 compute core per task
#SBATCH --time=24:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=script7-generation-mmn
#SBATCH --output=script7-v2.out

module purge
source ./venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

python ./scripts/script7b-deterministic-MMN.py