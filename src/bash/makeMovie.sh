#!/bin/bash

#movie.sh
#Create a movie from files in separate directories

#arg1 - sim
#arg2 - name (of directory and movie)

#e.g. ./movie.sh 50A vr

#Arguments
sim="test"
name="hist"

#:)
folder=movie_files_${name}
#cd /home/oge1/lammps/sapphire/analysis/results/$sim

#Create directories
echo "Creating directories"
mkdir $folder

#Copy images to one folder
echo "Copying images"
cp img/$name/*.png $folder

#Create movie
echo "Creating movie"
/home/oge1/software/bin/ffmpeg -framerate 10 -i ${folder}/%*.png -vcodec mpeg4 -b 800k -r 10  movies/${name}_${sim}.avi -y

#Delete temporary files & directory
echo "Deleting temporary files"
rm -rf $folder

#Done! Good job :D
echo "Done!"

