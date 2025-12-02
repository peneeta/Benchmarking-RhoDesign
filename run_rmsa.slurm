#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH -t 10:00:00
#SBATCH --gpus=v100-32:2
#SBATCH --mail-user=pawojcik@andrew.cmu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# load modules
module load anaconda3/2024.10-1

# nav to PDB folder
cd ${HOME}/gen_ai_project_dir/RhoDesign/src

# temp1 seqs
# python ${HOME}/Benchmarking-RhoDesign/generated_seqs