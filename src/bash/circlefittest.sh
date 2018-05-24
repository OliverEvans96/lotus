#!/bin/bash
#circlefittest.sh
#Compile and run circlefittest.cpp

d="/home/oge1/lammps/sapphire/analysis/root"

echo "Compiling..."
if g++ `root-config --glibs --cflags` -lMathMore -lMinuit $d/circlefittest.cpp $d/CircleFitClass.cpp -o $d/circlefittest.out
then
	echo "Running..."
	./circlefittest.out
	echo "Done!"
else
	echo "Compilation failed - terminating"
fi
