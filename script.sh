#!/bin/bash
max=17
rm -Rf DGBTRF.dat DGBTRI.dat DGBSV.dat
for ((i = 4, j = 1 ; j<=$max ; i*=2, j++)); do
    echo "$j / $max , arg = $i"
    output=$( make -B tpPoisson1D_direct; taskset -c 7 ./bin/tpPoisson1D_direct $i )
    echo $i $(echo "$output" | grep -i dgbtrf) >> DGBTRF.dat
    echo $i $(echo "$output" | grep -i tridiag) >> DGBTRI.dat
    echo $i $(echo "$output" | grep -i dgbsv) >> DGBSV.dat
done

gnuplot -c plot.gp DGBTRF
gnuplot -c plot.gp DGBTRI
gnuplot -c plot.gp DGBSV
gnuplot plotall.gp
