set terminal pdf 
set output "rm3.pdf"

set samples 100000
set title "Heave of Float object vs Time (Multibody)"

set xlabel "Time (s)"
set xrange [0 : 40]

plot "output.txt" using 1:2 with lines