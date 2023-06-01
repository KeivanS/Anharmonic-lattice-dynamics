set terminal postscript enhanced color
set output "dosplot"
set title "Density of states"
set xlabel "frequency(/cm)"
set key outside right
plot "./dos_tet.dat" u 2:4 w l lw 2 lt 1 noti
!ps2pdf dosplot
system("rm dosplot")
