set terminal pdf 
set output "oswec.pdf"

set samples 100000
set grid

# set lmargin at screen 0.125
set ylabel "Flap Rotation about y (degrees)"
set xrange [0:300]
plot "oswec_decay.txt" using 1:3 with lines title "most recent run"

# plot"oswec_all.txt" using 1:3 with lines title "all hydro forces",\
#      "oswec_no_radiation.txt" using 1:3 with lines title "no radiation damping force",\
#      "oswec_decay.txt" using 1:3 with lines title "most recent run"
    # "oswec_no_hydrostatic.txt" using 1:3 with lines title "no hydrostatic force", \
#  "oswec_none.txt" using 1:3 with lines title "no hydro forces", \
