set terminal pdf size 5in,4.5in
set output "rm3.pdf"

set samples 100000
set xrange [0:10]

set lmargin at screen 0.125
set multiplot layout 2,1
set ylabel "Float Heave"
set format x ''
plot "rm3_decay.txt" using 1:2 with lines title "Chrono", \
     "rm3_wec-sim_verification_data.txt" using 1:2 with lines title "WEC-Sim"

set ylabel "Plate Heave"
set xlabel "Time (s)"
unset format x 
plot "rm3_decay.txt" using 1:3 with lines title "Chrono", \
     "rm3_wec-sim_verification_data.txt" using 1:3 with lines title "WEC-Sim"