#!/bin/bash

cd ../results/50A

#Copy all quiver pics
mkdir -p allVr
cp atom*/img/vr/*.png allVr

#Create movie
/home/oge1/bin/ffmpeg -framerate 10 -i allVr/%*.png -vcodec mpeg4 -b 800k -r 30  allVr.avi -y
