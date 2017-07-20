#!/bin/bash

#Run quiverTest for 50A/atom1
#Execute from analysis/root

whichone="40A/atom2"

#Copy center of mass data
cp ../results/$whichone/center_of_mass.txt .

echo "Compiling . . ."
if ./compilePS.sh 
then
	echo "Compilation successful - Executing!"
	./polarScatter.out /home/oge1/lammps/sapphire/analysis/results/${whichone}/calculated.txt ${whichone} $1
else
	echo "Compilation failed - Terminating"
fi
