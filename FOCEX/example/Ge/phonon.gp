set terminal postscript enhanced color
set output "phononplot"
set title "Phonon dispersion"
set ylabel "frequency(/cm)"
set key outside right
plot "./mode_grun.dat" u 3:7 w p pt 7 ps 0.4 lc rgb "blue" noti
!ps2pdf phononplot
system("rm phononplot")
