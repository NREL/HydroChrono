set terminal pdf size 6in,4.5in
set output "deepcwind_decay.pdf"

set samples 100000
set grid
set lmargin at screen 0.125
set multiplot title "DeepCwind" # layout 2,1
set xlabel "Time (s)"
set xrange [0 : 1000]

# f(x) = c 
# fit f(x) 'deepcwind_LC44.txt' using 1:2 every ::1 via c 

# set title "Heave Decay (LC 4.4)"
# plot "NREL_LC_44.txt" every ::1 using 1:4 with lines linewidth 2.5 title "Compare", \
#      "deepcwind_LC44.txt" every ::1 using 1:2 with lines title "Chrono", \
#      f(x) with lines linewidth 0.7 title "equilibrium"

# f(x) = c 
# fit f(x) 'deepcwind_LC42.txt' using 1:2 every ::1 via c 

# set title "Surge Decay (LC 4.2)"
# plot "NREL_LC_42.txt" every ::1 using 1:2 with lines linewidth 2.5 title "Compare", \
#      "deepcwind_LC42.txt" every ::1 using 1:2 with lines title "Chrono", \
#      f(x) with lines linewidth 0.7 title "equilibrium"
     
f(x) = c 
fit f(x) 'deepcwind_LC46.txt' using 1:2 every ::1 via c 
set format y "%g°" 
set title "Pitch Decay (LC 4.6)"
plot "NREL_LC_46.txt" every ::1 using 1:6 with lines linewidth 2.5 title "Compare", \
     "deepcwind_LC46.txt" every ::1 using 1:($2 * 180 / 3.1415) with lines title "Chrono", \
     f(x) with lines linewidth 0.7 title "equilibrium"