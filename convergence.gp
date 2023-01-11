
set term png
set output "img/convergence.png"
set title "Convergence du residu par la méthode de Richardson"
set ylabel "Residu relatif"
set xlabel "Iterations"
set datafile separator " "
plot "data/RESVECalpha.dat" title "richardson_alpha" with lines, \
        "data/RESVECJac.dat" title "Jacobi" with lines, \
        "data/RESVECGS.dat" title "Gauss-Seidel" with lines
