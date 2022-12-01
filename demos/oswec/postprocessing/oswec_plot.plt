set terminal pdf 
set output "oswec.pdf"

set samples 100000
set grid

# set lmargin at screen 0.125
set ylabel "Flap Rotation about y (Radians)"
set xrange [0:500]
plot "oswec_decay.txt" using 1:2 with lines title "Chrono",\
     "wecsim_comp.txt" using 1:2 with lines title " WECSim"