set xlabel "Pair distance (Ang)"
set ylabel "Trace of FC2  (eV/Ang^2)"
set title "FC2 strangth (different color for different pair interactions)"
set yzeroaxis
set grid xtics
set grid ytics

set palette maxcolors 10
set palette defined   
set cblabel "Primitive atom"
set cbtics
#set colorbox user origin 0.92,0.05  # Places box at 95% right, 5% up from bottom
#set colorbox size 0.02,0.85     
set cbrange [1:5 ] # same range as palette definition

### plot [1:] 'trace_fc.dat' u 1:2 w p pt 1 ps 2.8 title 'Original fc2' ,\

 plot [1:] 'trace_fc.dat' u 4:5:($1+$2) w p pt 1 ps 1.8 palette z title 'Original' ,\
       'trace_fc2_lr.dat' u 4:5:($1+$2) w p pt 4 ps 1.8 palette z title 'fc2 LR' ,\
 0 w l lt 5  lw 2 notitle
 pause -1
     # 'trace_fc2_sr.dat' u 4:5 w p pt 9 ps 0.9 title 'fc2 SR'
#set arrow from x1, y1, z1 to x2, y2, z2 nohead lt 1
#set xrange [0:xmax]
#set yrange [0:ymax]
#set zrange [0:zmax]
#
# The file containing the arrow coordinates

reset

set style fill solid 1.0 border
set xlabel "X (Ang)"
set ylabel "Y (Ang)"
set zlabel "Z (Ang)"
file0 = 'atompos.xyz'
file1 = 'springs1.dat'
file2 = 'rgrid_raw.xyz'
file3 = 'supercell_1.xyz'
file4 = 'WSR_boundary.xyz'
file5 = 'trace_fc.dat'
# Plot the arrows
splot file0 u 2:3:4 with p pt 9 ps 0.7 title 'Atompos'   ,\
      file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
splot file1 u 1:2:3:4:5:6 with vectors nohead lt 1 title "1-j"  ,\
      file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file5 u 6:7:8:($9-$6):($10-$7):($11-$8) with vectors lt 3 title 'trace fc2' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
file1 = 'springs2.dat'
splot file1 u 1:2:3:4:5:6 with vectors nohead lt 1 title "2-j"  ,\
      file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
