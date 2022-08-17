set terminal pdf 
set output "sphere_decay.pdf"

set samples 100000
set title "Heave of Sphere vs Time (Multibody)"

set xlabel "Time (s)"
set xrange [0 : 40]

f(x) = c 
fit f(x) 'sphere_decay.txt' using 1:2 every ::1 via c 

plot "sphere_decay_comparison.txt" every ::1 using 1:($7-2) with lines linewidth 2.5 title "InWave+H (Lin) shifted down 2", \
     "sphere_decay.txt" every ::1 using 1:2 with lines title "Chrono", \
     f(x) with lines linewidth 0.7 title "Approximate new equilibrium"