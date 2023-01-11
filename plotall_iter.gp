
set term png
set output "ALL_iter.png"
set title "LAPACK functions"
set ylabel "Time in seconds"
set xlabel "Matrix size"
set datafile separator " "
plot "RICHARDSON.dat" using 1:5 with lines, \
         "JACOBI.dat" using 1:5 with lines, \
         "GAUSS.dat" using 1:5 with lines
