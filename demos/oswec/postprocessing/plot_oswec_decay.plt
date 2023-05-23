set terminal pdf 
set output "oswec_decay.pdf"

set samples 100000
set grid
set title "Oswec Pitch Decay Test"
set format y "%g°"
set yrange [-11:11]

# set lmargin at screen 0.125
set ylabel "Flap Pitch"
set xlabel "Time (s)"
set xrange [0:400]
plot "oswec_decay.txt" using 1:($2/3.14159*180.0) with lines title "Chrono",\
     "wecsim_oswec_decay.txt" using 1:($2/3.14159*180.0) with lines title " WECSim"