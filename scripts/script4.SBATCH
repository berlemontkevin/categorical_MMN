#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=10:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=script3

module purge
module load julia/1.5.3

~/julia/my-julia script4-strength-to-sst.jl