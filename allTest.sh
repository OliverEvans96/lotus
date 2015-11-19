#!/bin/bash

#Run quiverTest for 50A/atom1
#Execute from analysis/root
echo "Compiling . . ."
g++ `root-config --glibs --cflags` polarScatter.cpp CircleFitClass.cpp Quiver/Quiver.cpp -o polarScatter.out

echo "Executing!"
../exec/parallel.sh ../exec/analyze.sh Quiver50A_atom{} n004 /home/oge1/lammps/sapphire/analysis/results/50A/atom{}/calculated.txt 50A/atom{} polarScatter 1 12

cd ../results/50A

#Copy all quiver pics
mkdir -p allQuiver
cp atom*/quiver/*.png allQuiver

#Create movie
/home/oge1/bin/ffmpeg -framerate 10 -i allQuiver/%*.png -vcodec mpeg4 -b 800k -r 30  allQuiver.avi -y

