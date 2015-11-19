ffmpeg -framerate 10 -i img/step%07d.bmp -c:v libx264 -r 30 -pix_fmt yuv420p movie.mov
