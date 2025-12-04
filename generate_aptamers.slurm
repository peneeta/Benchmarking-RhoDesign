#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH -t 10:00:00
#SBATCH --gpus=v100-32:2
#SBATCH --mail-user=pawojcik@andrew.cmu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# load modules
module load anaconda3/2024.10-1
module load cuda/12.4.0

# activate python env
conda activate RhoDesign

# nav to PDB folder
cd ${HOME}/gen_ai_project_dir/RhoDesign/src

# apt 1
python inference_without2d.py -pdb ${HOME}/Benchmarking-RhoDesign/aptamers/2AU4.pdb -save ${HOME}/Benchmarking-RhoDesign/generated_seqs/2AU4 -temp 1

# apt 2
python inference_without2d.py -pdb ${HOME}/Benchmarking-RhoDesign/aptamers/5KPY.pdb -save ${HOME}/Benchmarking-RhoDesign/generated_seqs/5KPY -temp 1

# apt 3
if [ ! -d ${HOME}/Benchmarking-RhoDesign/generated_seqs/6E8T ]; then
  mkdir -p ${HOME}/Benchmarking-RhoDesign/generated_seqs/6E8T;
fi


python inference_without2d.py -pdb ${HOME}/Benchmarking-RhoDesign/aptamers/6E8T_mango.pdb -save ${HOME}/Benchmarking-RhoDesign/generated_seqs/6E8T -temp 1
