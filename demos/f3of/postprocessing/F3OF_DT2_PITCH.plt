set terminal pdf
set output "F3OF_DT2.pdf"

set samples 100000
set grid
set lmargin at screen 0.125
set multiplot title "F3OF DT2 Pitch"
set xrange [0:300]
set format y "%g°"

set ylabel "Pitch"
# set key at screen 1, 0.55
set key rmargin
set title " "
plot "reference_data/INW/DT1-3/INW_F3OF_DT2_PITCH.dat" using 1:2 with lines lt rgb "yellow4" title "INW", \
     "reference_data/PDS/DT1-3/PDS_F3OF_DT2_PITCH.dat" using 1:2 with lines lt rgb "light-red" title "PDS", \
     "reference_data/WDN/DT1-3/WDN_F3OF_DT2_PITCH.dat" using 1:2 with lines title "WDN", \
     "reference_data/WSM/DT1-3/WSM_F3OF_DT2_PITCH.dat" using 1:2 with lines title "WSM", \
     "CHRONO_F3OF_DT2_PITCH.txt" using 1:($3/3.14159*180.0) with lines lt rgb "dark-violet" title "Chrono"