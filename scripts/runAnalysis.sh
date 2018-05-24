#!/bin/bash

#Run quiverTest for 50A/atom1
#Execute from analysis/root
echo "Compiling . . ."
g++ `root-config --glibs --cflags` polarScatter.cpp CircleFitClass.cpp Quiver.cpp -o polarScatter.out

echo "Executing!"
../exec/analyze.sh /home/oge1/lammps/sapphire/analysis/results/50A/atom1/calculated.txt 50A/atom1 polarScatter

