#!/usr/bin/env gnuplot -c

# 1. Handle arguments
Th = (ARGC > 0) ? ARG1 + 0.0 : 1.0
N  = (ARGC > 1) ? ARG2 : "40"

# 2. Run simulation
cmd = sprintf("./bandfs_reaction_diffusion %s %f > result.txt", N, Th)
system(cmd)

# 3. Aesthetics
set title sprintf("Steady-State Reaction-Diffusion ({/Symbol f} = %.2f)", Th)
set xlabel "Dimensionless Position (x)"
set ylabel "Concentration u(x)"
set grid
set key left top
set yrange [0:1.05]

# Tell gnuplot to ignore lines starting with 'x' (the header)
set datafile commentschars "x"

# 4. Analytical Solution
f(x) = cosh(Th * x) / cosh(Th)

# 5. Plot
plot "result.txt" u 1:2 w lp pt 7 title "Numerical (FD)", \
     f(x) title "Analytical"