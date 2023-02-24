set terminal pdf 
set output "sphere_decay.pdf"
set termoption dashed

set samples 100000
set title "Sphere Heave Decay Test"

set grid
set xlabel "Time (s)"
set ylabel "Heave (m)"
set xrange [0 : 40]

f(x) = c 
fit f(x) 'sphere_decay.txt' using 1:2 every ::1 via c 

plot "sphere_decay.txt" every ::1 using 1:2 with lines linewidth 2.5 title "Chrono", \
     "sphere_decay_comparison.txt" every ::1 using 1:($7-2) with lines dashtype 4 linetype 5 title "InWave+H (Lin)", \
     f(x) with lines linewidth 0.7 title "Calculated Equilibrium"