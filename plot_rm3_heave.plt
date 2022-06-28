set terminal pdf 
set output "rm3_heave.pdf"

set samples 100000
set title "Heave of Float object vs Time"

set xlabel "Time (s)"
set xrange [1 : 40]

plot "output.txt" using 1:2 with lines