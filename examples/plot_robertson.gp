
set terminal pngcairo enhanced font ",12"
set output "robertson_rate.png"

#set logscale x

set title "Robertson rate equations"

set xlabel "Time"
set ylabel "Concentration"

plot "r1.txt" u 1:2 w lp lw 2 title "c_A",\
     "r1.txt" u 1:($3*1.E4) w lp lw 2 title "10^4 × c_B",\
     "r1.txt" u 1:4 w lp lw 2 title "c_C",\
