set terminal postscript enhanced color
set output "phononplot"
set title "Phonon dispersion with varying neighbors}"
set ylabel "frequency(/cm)}"
set key outside right
plot "./Ge-11105/mode_grun.dat" u 3:7 w p pt 6 ps 0.8 ti "1105","./Ge-11106/mode_grun.dat" u 3:7 w p pt 6 ps 0.8 lc rgb "red" ti "1106","./Ge-11107/mode_grun.dat" u 3:7 w p pt 6 ps 0.8 lc rgb "green" ti "1107","./Ge-11108/mode_grun.dat" u 3:7 w p pt 6 ps 0.8 lc rgb "yellow" ti "1108","./Ge-11109/mode_grun.dat" u 3:7 w p pt 6 ps 0.8 lc rgb "brown" ti "1109","../../THERMACOND/old_log/mode_grun.dat" u 3:7 w p pt 6 ps 0.8 lc rgb "black" ti "Ref"
!ps2pdf phononplot
