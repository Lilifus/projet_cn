set term png
set output sprintf("./%s.png",ARG1)
set title "LAPACK functions"
set ylabel "Time in seconds"
set xlabel "Matrix size"
set datafile separator " "
plot sprintf("./%s.dat",ARG1) using 1:3 with lines
