set term png
set output sprintf("img/%s.png",ARG1)
set title "LAPACK functions"
set ylabel "Time in seconds"
set xlabel "Matrix size"
set datafile separator " "
plot sprintf("data/%s.dat",ARG1) using 1:3 title sprintf("%s",ARG1) with lines
