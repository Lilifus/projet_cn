
set term png
set output "img/ALL_direct.png"
set title "LAPACK functions"
set ylabel "Time in seconds"
set xlabel "Matrix size"
set datafile separator " "
plot "data/DGBTRF.dat" using 1:3 title "dgbtrf" with lines, \
         "data/DGBTRI.dat" using 1:3 title "dgbtridiag" with lines, \
         "data/DGBSV.dat" using 1:3 title "dgbtrsv" with lines
