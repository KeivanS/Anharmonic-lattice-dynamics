set xlabel "Pair distance (Ang)"
set ylabel "Trace of FC2 (eV/Ang^2)"
set yzeroaxis
set grid xtics
 plot [1:] 'trace_fc.dat' w p notitle
 pause -1
#set arrow from x1, y1, z1 to x2, y2, z2 nohead lt 1
#set xrange [0:xmax]
#set yrange [0:ymax]
#set zrange [0:zmax]
#
# The file containing the arrow coordinates
set style fill solid 1.0 border
unset xlabel
unset ylabel
file0 = 'atompos.xyz'
file1 = 'springs1.dat'
file2 = 'rgrid_raw.xyz'
file3 = 'supercell_1.xyz'
file4 = 'WSR_boundary.xyz'
# Plot the arrows
splot file0 u 2:3:4 with p pt 9 ps 0.7 title 'Atompos'   ,\
      file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
splot file1 u 1:2:3:4:5:6 with vectors nohead lt 1 title "1-j"  ,\
      file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
file1 = 'springs2.dat'
splot file1 u 1:2:3:4:5:6 with vectors nohead lt 1 title "2-j"  ,\
      file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
file1 = 'springs3.dat'
splot file1 u 1:2:3:4:5:6 with vectors nohead lt 1 title "O-j"  ,\
      file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
