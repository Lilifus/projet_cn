
set term png size 800,600
set output "img/conv_richardson.png"
set title "Convergence du residu par la méthode de Richardson"
set ylabel "Residu relatif"
set xlabel "Iterations"
set datafile separator " "
plot "data/RESVECalpha.dat" title "Richardson_a_l_p_h_a" with lines

set term png size 800,600
set output "img/conv_jacobi.png"
set title "Convergence du residu par la méthode de Jacobi"
set ylabel "Residu relatif"
set xlabel "Iterations"
set datafile separator " "
plot "data/RESVECJac.dat" title "Jacobi" with lines

set term png size 800,600
set output "img/conv_gauss.png"
set title "Convergence du residu par la méthode de Gauss-Seidel"
set ylabel "Residu relatif"
set xlabel "Iterations"
set datafile separator " "
plot "data/RESVECGS.dat" title "Gauss-Seidel" with lines
