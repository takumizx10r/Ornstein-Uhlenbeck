#!/bin/sh
RAMDA=0.00
#for RAMDA in 0.0 0.01 0.02 0.05 0.1 0.2 0.5 1.0; do
for RAMDA in 0.0 0.01 1.0; do
echo "${RAMDA}"
 g++ -o run main.cpp -lm -Wall -DRAMDA=${RAMDA} 
 ./run
 gnuplot gnuplot/animation.gp
  ffmpeg -y -framerate 100 -i gnuplot/image/image%05d.png\
 -pix_fmt yuv420p "gnuplot/movies/movie_${RAMDA}.mp4"
# gnuplot gnuplot/draw-msd.gp
# gnuplot gnuplot/draw-diffusionconst.gp
 gnuplot -e "RAMDA=${RAMDA}" gnuplot/draw-msd.gp
 gnuplot -e "RAMDA=${RAMDA}" gnuplot/draw-diffusionconst.gp
done

#  gnuplot -e "K=${K}" gnuplot/animationwithhistgram.gp
#  
#gnuplot gnuplot/temp.gp
