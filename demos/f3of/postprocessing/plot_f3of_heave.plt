set terminal pdf size 5in,4.5in
set output "f3of_decay.pdf"

set samples 100000
set grid

set lmargin at screen 0.125
set multiplot layout 2,1
set title "F3OF Decay test DT3"
# set ylabel "Base Heave"
set xrange [0:300]
# plot "f3of_decay.txt" using 1:2 with lines title "Chrono"

set ylabel "Flap Pitch (radians)"
set xlabel "Time (s)"
set ylabel "Fore Flap"
plot "f3of_decay.txt" using 1:3 with lines title "HydroChrono"
set ylabel "Aft Flap"
plot "f3of_decay.txt" using 1:4 with lines title "HydroChrono"