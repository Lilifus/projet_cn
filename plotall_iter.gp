
set term png
set output "img/ALL_iter.png"
set title "LAPACK functions"
set ylabel "Time in seconds"
set xlabel "Matrix size"
set datafile separator " "
plot "data/RICHARDSON.dat" using 1:5 title "richardson_alpha" with lines, \
         "data/JACOBI.dat" using 1:5 title "jacobi" with lines, \
         "data/GAUSS.dat" using 1:5 title "gauss-seidel" with lines
