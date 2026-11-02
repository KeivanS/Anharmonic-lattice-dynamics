set xlabel "Pair distance (Ang)"
set ylabel "Trace of FC2  (eV/Ang^2)"
set title "FC2 strength (Contour 2 for 1-1 type, 3 for 1-2 type and 4 for 2-2 type pair interactions)"
set yzeroaxis
set grid xtics
set grid ytics
set key bottom right
#set palette maxcolors 10
set palette defined   
set cblabel "Primitive atom"
set cbtics
#set colorbox user origin 0.92,0.05  # Places box at 95% right, 5% up from bottom
#set colorbox size 0.02,0.85     
set cbrange [1:5 ] # same range as palette definition

### plot [1:] 'trace_fc.dat' u 1:2 w p pt 1 ps 2.8 title 'Original fc2' ,\

 plot [1:] 'trace_fc.dat' u 4:5:($1+$2) w p pt 1 ps 1.8 palette z title 'Original' ,\
       'trace_phi_lr.dat' u 4:5:($1+$2) w p pt 4 ps 1.8 palette z title 'fc2 on the extended grid' ,\
 0 w l lt 5  lw 2 notitle
 pause -1
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
set view equal xyz
set ticslevel 0
file0 = 'atompos.xyz'
file2 = 'rgrid_raw.xyz'
file3 = 'supercell_1.xyz'
file4 = 'WSR_boundary.xyz'
splot file0 u 2:3:4 with p pt 9 ps 0.7 title 'Atompos'   ,\
      file2 u 2:3:4 with p pt 7 ps 0.6 title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.5 title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.1 title 'WScell of SC'
pause -1
# one plot per primitive atom i: its FC2 bonds (springs<i>.dat: r_i - r_j) drawn from atom i to each
# neighbour j, and those neighbours
nat = int(system("ls springs[0-9]*.dat 2>/dev/null | wc -l"))
do for [i=1:nat] {
    fs = sprintf('springs%d.dat', i)
    pos = system(sprintf("awk '$6==\"%d(\" && $7==0 && $8==0 && $9==\"0)\" {print $2, $3, $4; exit}' %s", i, file0))
    xi = real(word(pos, 1)); yi = real(word(pos, 2)); zi = real(word(pos, 3))
    set title sprintf("FC2 neighbours of atom %d", i)
    splot fs u (xi):(yi):(zi):(-$4):(-$5):(-$6) with vectors nohead lt 1 title sprintf("%d-j", i) ,\
          fs u (xi-$4):(yi-$5):(zi-$6) with p pt 7 ps 0.5 lc rgb "black" title 'neighbours' ,\
          '+' u (xi):(yi):(zi) every ::0::0 with p pt 7 ps 1.5 lc rgb "red" title sprintf('atom %d', i) ,\
          file4 u 2:3:4 with p pt 2 ps 0.1 lc rgb "gray" title 'WScell of SC'
    pause -1
}
