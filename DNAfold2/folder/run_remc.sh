#!/bin/bash
#SBATCH -n 10
rm -rf model
mkdir model
bash initial.sh
cp DNAfold2_remc.c config1.dat model
cd model
gcc -O3 -Wall -fopenmp DNAfold2_remc.c -o DNAfold2_remc -lm
./DNAfold2_remc
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
g++ tm.c
./a.out
