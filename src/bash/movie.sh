ffmpeg -framerate 10 -i img/$1/%*.png -vcodec mpeg4 -b 800k -r 10 img/$1.avi
