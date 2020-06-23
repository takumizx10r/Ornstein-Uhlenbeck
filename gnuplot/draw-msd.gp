reset
set terminal pngcairo enhanced size 1200,1200 font "Arial, 28"
D0=0.0100000
inputfile = sprintf("dataset/data%.4f.dat",RAMDA)
outputfile =sprintf("gnuplot/msd/msd%.4f.png",RAMDA)
#############################################################
f(x) =  4.0*D0*x**a
a=0.1
TITLE = sprintf("RAMDA=%.4f",RAMDA)
set output outputfile
unset xrange
set yrange[0.010:100]
set log
plot inputfile u 1:2 pt 7 title TITLE, inputfile u 1:(4*$4*$1) with lines title "pure-"
#fit f(x) inputfile u 1:2 via a
#plot inputfile u 1:2 pt 7, inputfile u 1:(4*$4*$1) with lines, f(x) with lines lc "black"

##############################################################
#set print savefile append
#print K,a,b
