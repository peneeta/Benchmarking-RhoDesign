#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH -t 10:00:00
#SBATCH --gpus=v100-32:2
#SBATCH --mail-user=pawojcik@andrew.cmu.edu
#SBATCH --mail-type=BEGIN,END,FAIL


# TODO add code to make dirs

# load modules
module load anaconda3/2024.10-1
module load cuda/12.4.0

# activate python env
conda activate RhoDesign

# nav to PDB folder
cd ${HOME}/gen_ai_project_dir/RhoDesign/src

# apt 1
if [ ! -d ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/2AU4 ]; then
  mkdir -p ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/2AU4;
fi

python inference_without2d.py -pdb ${HOME}/Benchmarking-RhoDesign/aptamers/2AU4.pdb -save ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/2AU4 -temp 1e-5

# apt 2
if [ ! -d ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/5KPY ]; then
  mkdir -p ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/5KPY;
fi

python inference_without2d.py -pdb ${HOME}/Benchmarking-RhoDesign/aptamers/5KPY.pdb -save ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/5KPY -temp 1e-5

# apt 3
if [ ! -d ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/6E8T ]; then
  mkdir -p ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/6E8T;
fi

python inference_without2d.py -pdb ${HOME}/Benchmarking-RhoDesign/aptamers/6E8T_mango.pdb -save ${HOME}/Benchmarking-RhoDesign/generated_seqs_lowtemp/6E8T -temp 1e-5

