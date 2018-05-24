#!/bin/bash
d="/home/oge1/lammps/sapphire/analysis/root"
g++ `root-config --glibs --cflags` -lMathMore $d/NumericalMinimization.cpp -o $d/NumericalMinimization.out
