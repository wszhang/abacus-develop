#!/bin/sh
module purge && module load anaconda3_nompi gcc/9.2.0 2>&1
module list 2>&1;
source activate pytorch110
pwd;
ls;
python3 ../SIAB.py ORBITAL_INPUT_DZP
