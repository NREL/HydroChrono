set terminal pdf size 6in,4.5in
set output "F3OF_DT3.pdf"

set samples 100000
set grid
set lmargin at screen 0.125
set multiplot title "F3OF DT3 Decay Test"
set xrange [0:300]
set format y "%g°"

set margins 10,14,0,3
set format x ""
set ylabel "Fore Flap"
set key at screen 1, 0.55
set size 1, .475
set origin 0, .525
set title " "
plot "CHRONO_F3OF_DT3_PITCH.txt" using 1:($4/3.14159*180.0) with lines title "Chrono", \
     "reference_data/INW/DT1-3/INW_F3OF_DT3_FLAP1.dat" using 1:2 with lines title "INW", \
     "reference_data/PDS/DT1-3/PDS_F3OF_DT3_FLAP1.dat" using 1:2 with lines title "PDS", \
     "reference_data/WDN/DT1-3/WDN_F3OF_DT3_FLAP1.dat" using 1:2 with lines title "WDN", \
     "reference_data/WSM/DT1-3/WSM_F3OF_DT3_FLAP1.dat" using 1:2 with lines title "WSM"

set margins 10,14,2,2
set format x "%g"
set ylabel "Aft Flap"
unset key
set size 1, .475
set origin 0, .05
unset title
set xlabel "Time (s)"
plot "CHRONO_F3OF_DT3_PITCH.txt" using 1:($5/3.14159*180.0) with lines title "Chrono", \
     "reference_data/INW/DT1-3/INW_F3OF_DT3_FLAP2.dat" using 1:2 with lines title "INW", \
     "reference_data/PDS/DT1-3/PDS_F3OF_DT3_FLAP2.dat" using 1:2 with lines title "PDS", \
     "reference_data/WDN/DT1-3/WDN_F3OF_DT3_FLAP2.dat" using 1:2 with lines title "WDN", \
     "reference_data/WSM/DT1-3/WSM_F3OF_DT3_FLAP2.dat" using 1:2 with lines title "WSM"