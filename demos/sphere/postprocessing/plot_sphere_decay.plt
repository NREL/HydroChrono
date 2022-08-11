set terminal pdf 
set output "sphere_decay.pdf"

set samples 100000
set title "Heave of Float object vs Time (Multibody)"

set xlabel "Time (s)"
set xrange [0 : 40]

plot "sphere_decay.txt" using 1:2 with lines title "Chrono", \
     "WECSimSphereDecayCompare.txt" using 1:2 with lines title "WEC-Sim"