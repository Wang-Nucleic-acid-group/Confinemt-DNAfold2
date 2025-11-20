#!/bin/bash
#SBATCH -n 1
rm -rf model
mkdir model
bash initial.sh
cp DNAfold2_sa.c config1.dat model
cd model
gcc -O3 -Wall DNAfold2_sa.c -o DNAfold2_sa -lm
./DNAfold2_sa
cd ..
cp t1.c model
cd model
g++ t1.c
./a.out
cd ..
bash scoring.sh
bash secondary.sh
bash optimize.sh
bash re_scoring.sh
bash rebuild.sh
bash tm.sh
