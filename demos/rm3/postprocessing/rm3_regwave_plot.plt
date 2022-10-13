set terminal pdf
set output "reg_waves.pdf"

set samples 100000
set grid

# set lmargin at screen 0.125
# set multiplot layout 2,1
set ylabel "Float Heave"
set xrange [0:40]
plot "regwave_10.txt" using 1:2 with lines title "Chrono"

# set ylabel "Plate Heave"
# set xlabel "Time (s)"
# set xrange [0:300]
# plot "rm3_decay.txt" using 1:3 with lines title "Chrono", \
#      "rm3_WECSim_decay.txt" using 1:3 with lines title "WEC-Sim"