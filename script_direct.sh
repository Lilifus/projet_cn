#!/bin/bash
max=25
rm -Rf data/DGBTRF.dat data/DGBTRI.dat data/DGBSV.dat
for ((i = 4, j = 1 ; j<=$max ; i*=2, j++)); do
    echo "$j / $max , arg = $i"
    output=$( make -B tpPoisson1D_direct; taskset -c 1 ./bin/tpPoisson1D_direct $i )
    echo $i $(echo "$output" | grep -i dgbtrf) >> data/DGBTRF.dat
    echo $i $(echo "$output" | grep -i tridiag) >> data/DGBTRI.dat
    echo $i $(echo "$output" | grep -i dgbsv) >> data/DGBSV.dat
done

mkdir -p img
gnuplot -c plot_direct.gp DGBTRF
gnuplot -c plot_direct.gp DGBTRI
gnuplot -c plot_direct.gp DGBSV
gnuplot plotall_direct.gp
