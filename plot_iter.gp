set term png size 800,600
set output sprintf("img/%s.png",ARG1)
set title "LAPACK functions"
set ylabel "Time in seconds"
set xlabel "Matrix size"
set datafile separator " "
plot sprintf("data/%s.dat",ARG1) using 1:5 title sprintf("%s",ARG1) with lines
