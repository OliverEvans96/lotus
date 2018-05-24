#!/bin/bash
d="/home/oge1/lammps/sapphire/analysis/root"
g++ `root-config --glibs --cflags` -lMathMore -lMinuit $d/cyl_analysis.cpp $d/CircleFitClass.cpp $d/Quiver/Quiver.cpp -o $d/cyl_analysis.out
