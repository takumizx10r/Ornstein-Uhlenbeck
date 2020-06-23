set terminal png enhanced size 1200, 1200 font "Arial, 18"
ROI=1.0
MAXITER=1000
# Format of the axis
if(exist("n") == 0) n = 0
set xrange [-ROI:ROI]
set yrange [-ROI:ROI]
set xlabel "x"
set ylabel "y"
set format x "%.2f"
set format y "%.2f"
#set key font "Arial, 18"
#set key spacing 1.5 
set size square 
inputFile  = sprintf("result/result%05d.dat",n)
#inputFileforwall  = sprintf("result/wall%05d.dat",n)
outputFile = sprintf("gnuplot/image/image%05d.png",n)
set output outputFile
plot inputFile with points pt 7 ps 0.5 lc rgb "#ff00ff" notitle

clear

n=n+1
if(n <= MAXITER) reread
