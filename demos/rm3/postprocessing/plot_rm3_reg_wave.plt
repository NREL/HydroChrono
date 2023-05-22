set terminal pdf size 5in, 5in
set output "reg_waves.pdf"
set encoding utf8

set samples 100000
set grid

set xlabel "Time (s)"
# set lmargin at screen 0.125
set multiplot layout 2,1 title "RM3 Regular Waves: A=0.022, ω=2.10"
set ylabel "Float Heave (m)"
set xrange [0:40]
plot "rm3_reg_waves.txt" using 1:2 with lines title "Chrono"
# , \
# "rm3_wec-sim_verification_data.txt" using 1:2 with lines title "WEC-Sim"

set ylabel "Drift (m)"
set xrange [0:40]
plot "rm3_reg_waves.txt" using 1:4 with lines title ""

# set ylabel "Plate Heave (m)"
# set xrange [0:40]
# plot "rm3_reg_waves.txt" using 1:3 with lines title "Chrono"
# , \
#      "rm3_wec-sim_verification_data.txt" using 1:3 with lines title "WEC-Sim"