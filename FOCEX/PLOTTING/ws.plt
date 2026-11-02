set yzeroaxis
set grid xtics
file2 = 'rgrid_raw.xyz'
file4 = 'WSR_boundary.xyz'
# Plot the arrows
splot file2 u 2:3:4 with p pt 7 ps 1.5 title 'R-vectors' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of SC'
pause -1
file2 = 'ggrid_raw.xyz'
file4 = 'WSG_boundary.xyz'
splot file2 u 2:3:4 with p pt 7 ps 1.5 title 'G-vectors' ,\
      file4 u 2:3:4 with p pt 2 ps 0.2 title 'WScell of BZ'
pause -1
