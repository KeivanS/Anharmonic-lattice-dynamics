set terminal postscript enhanced color
set output "dosplot"
set title "DOS with varying neighbors"
set xlabel "frequency(/cm)"
set key outside right
plot "./Ge-11105/dos_tet.dat" u 2:4 w l lw 2 lt 1 ti "1105","./Ge-11106/dos_tet.dat" u 2:4 w l lw 2 lt 1 lc rgb "red" ti "1106","./Ge-11107/dos_tet.dat" u 2:4 w l lw 2 lt 1 lc rgb "green" ti "1107","./Ge-11108/dos_tet.dat" u 2:4 w l lw 2 lt 1 lc rgb "yellow" ti "1108","./Ge-11109/dos_tet.dat" u 2:4 w l lw 2 lt 1 lc rgb "brown" ti "1109"
!ps2pdf dosplot
