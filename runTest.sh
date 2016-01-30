#!/bin/bash

#Run quiverTest for 50A/atom1
#Execute from analysis/root
echo "Compiling . . ."
./compilePS.sh

echo "Executing!"
./polarScatter.out /home/oge1/lammps/sapphire/analysis/results/40A/atom2/calculated.txt 40A/atom2 $1

#cd ../results/50A

#Copy all quiver pics
#mkdir -p allQuiver
#cp atom*/quiver/*.png allQuiver

#Create movie
#/home/oge1/bin/ffmpeg -framerate 10 -i allQuiver/%*.png -vcodec mpeg4 -b 800k -r 30  allQuiver.avi -y

