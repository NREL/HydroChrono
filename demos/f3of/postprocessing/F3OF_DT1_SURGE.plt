set terminal pdf
set output "F3OF_DT1.pdf"

set samples 100000
set grid
set lmargin at screen 0.125
set multiplot title "F3OF DT1 Surge Decay Test"
set xrange [0:300]
# set format y "%g°"

# set margins 10,14,0,3
# set format x ""
# unset xlabel
set ylabel "Surge (m)"
set xlabel "Time (s)"
# set key at screen 1, 0.55
set key rmargin
# set size 1, .5
# set origin 0, .5
set title " "
plot "CHRONO_F3OF_DT1_SURGE.txt" using 1:2 with lines title "Chrono", \
     "reference_data/INW/DT1-3/INW_F3OF_DT1_SURGE.dat" using 1:2 with lines title "INW", \
     "reference_data/PDS/DT1-3/PDS_F3OF_DT1_SURGE.dat" using 1:2 with lines title "PDS", \
     "reference_data/WDN/DT1-3/WDN_F3OF_DT1_SURGE.dat" using 1:2 with lines title "WDN", \
     "reference_data/WSM/DT1-3/WSM_F3OF_DT1_SURGE.dat" using 1:2 with lines title "WSM"

# set margins 10,14,2,2
# set format x "%g"
# set xlabel "Time (s)"
# set ylabel "Aft Flap"
# unset key
# set size 1, .5
# set origin 0, 0
# unset title
# plot "CHRONO_F3OF_DT1_PITCH.txt" using 1:($4/3.14159*180.0) with lines title "Chrono", \
#      "reference_data/INW/DT1-3/INW_F3OF_DT1_FLAP2.dat" using 1:2 with lines title "INW", \
#      "reference_data/PDS/DT1-3/PDS_F3OF_DT1_FLAP2.dat" using 1:2 with lines title "PDS", \
#      "reference_data/WDN/DT1-3/WDN_F3OF_DT1_FLAP2.dat" using 1:2 with lines title "WDN", \
#      "reference_data/WSM/DT1-3/WSM_F3OF_DT1_FLAP2.dat" using 1:2 with lines title "WSM"