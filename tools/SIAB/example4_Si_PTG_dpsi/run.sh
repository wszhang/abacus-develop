#!/bin/sh
module purge && module load anaconda3_nompi gcc/9.2.0 elpa/2021.05.002/intelmpi2018 intelmpi/2018.update4 2>&1
module list 2>&1;
source activate pytorch110
pwd;
ls;
python3 ../SIAB.py ORBITAL_INPUT_DZP

