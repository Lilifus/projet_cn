
set term png
set output "ALL.png"
set title "LAPACK functions"
set ylabel "Time in seconds"
set xlabel "Matrix size"
set datafile separator " "
plot "DGBTRF.dat" using 1:3 with lines, \
         "DGBTRI.dat" using 1:3 with lines, \
         "DGBSV.dat" using 1:3 with lines
