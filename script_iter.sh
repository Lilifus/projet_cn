#!/bin/bash
max=14
rm -Rf RICHARDSON.dat JACOBI.dat GAUSS.dat
for ((i = 4, j = 1 ; j<=$max ; i*=2, j++)); do
    echo "$j / $max , arg = $i"
    output=$( make -B tpPoisson1D_iter; taskset -c 1 ./bin/tpPoisson1D_iter $i )
    echo $i $(echo "$output" | grep -i richardson_alpha) >> RICHARDSON.dat
    echo $i $(echo "$output" | grep -i jacobi) >> JACOBI.dat
    echo $i $(echo "$output" | grep -i gauss) >> GAUSS.dat
done

gnuplot -c plot.gp RICHARDSON
gnuplot -c plot.gp JACOBI
gnuplot -c plot.gp GAUSS
gnuplot plotall_iter.gp
