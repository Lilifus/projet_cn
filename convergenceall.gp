
set term png size 800,600
set output "img/convergence.png"
set title "Convergence du residu par la méthode de Richardson"
set ylabel "Residu relatif"
set xlabel "Iterations"
set datafile separator " "
plot "data/RESVECalpha.dat" title "Richardson_a_l_p_h_a" with lines, \
        "data/RESVECJac.dat" title "Jacobi" with lines, \
        "data/RESVECGS.dat" title "Gauss-Seidel" with lines
