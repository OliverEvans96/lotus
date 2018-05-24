#!/bin/bash

#Run quiverTest for 50A/atom1
#Execute from analysis/root

whichone="Bob/Sub951By50/Cyl30A/atom1"

#Copy center of mass data
cp ../results/$whichone/center_of_mass.txt .

echo "Compiling . . ."
if ./compileCYL.sh 
then
	echo "Compilation successful - Executing!"
	./cyl_analysis.out /home/oge1/lammps/sapphire/analysis/results/${whichone}/calculated.txt ${whichone} $1
else
	echo "Compilation failed - Terminating"
fi
