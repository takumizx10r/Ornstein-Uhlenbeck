reset
set terminal pngcairo enhanced size 1200,1200 font "Arial, 28"
D0=0.0100000
inputfile  =  sprintf("dataset/data%.4f.dat",RAMDA)
outputfile =  sprintf("gnuplot/diffusion-constant/d-c%.4f.png",RAMDA)
#############################################################
set output outputfile
set yrange [0:1.200]
set log x
TITLE = sprintf("RAMDA=%.4f",RAMDA)
plot inputfile u 1:($3/$4) with points title TITLE
