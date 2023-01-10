
set title "Convergence du residu par la méthode de Richardson"
set ylabel "Residu relatif"
set xlabel "Iterations"
set datafile separator " "
plot "RESVECalpha.dat" with lines, \
        "RESVECJac.dat" with lines, \
        "RESVECGS.dat" with lines
pause -1 "Press ENTER to continue"
