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

# trace_fc.dat: block (gnuplot index) 0 = FC2 traces, 1 = FC3 invariants, 2 = FC4 invariants
 plot [1:] 'trace_fc.dat' index 0 u 4:5:($1+$2) w p pt 1 ps 1.8 palette z title 'Original' ,\
       'trace_phi_lr.dat' u 4:5:($1+$2) w p pt 4 ps 1.8 palette z title 'fc2 on the extended grid' ,\
 0 w l lt 5  lw 2 notitle
 pause -1

# FC3 and FC4 rotational invariants per triplet/quadruplet versus its largest pair distance (2026-11-06,
# user), coloured by the number of distinct atoms: 2 = pair terms (i,i,j), 3 and 4 = many-body terms.
# Columns: 5 dmax, 6 Frobenius norm, 7 isotropic part (FC4), 8 number of distinct atoms. Plotted only
# when the block has data lines (the rank was included).
nblk(b) = int(system(sprintf("awk 'NF==0{e++;next} {if(e>=2)k++; e=0} k==%d && $1!~/^#/{n++} END{print n+0}' trace_fc.dat", b)))
n3 = nblk(1)
n4 = nblk(2)
# one graph, separate from FC2: ||FC3|| (circles) and ||FC4|| (triangles); then the signed FC4 isotropic part
if (n3 + n4 > 0) {
    reset
    set grid
    set logscale y
    set format y "10^{%T}"
    set xlabel "Largest pair distance in the triplet / quadruplet (Ang)"
    set ylabel "||FC3|| (eV/Ang^3),  ||FC4|| (eV/Ang^4)"
    set title "FC3 and FC4 strength: sqrt(sum of the squared Cartesian components) per cluster"
    set key top right
    set palette maxcolors 4
    set palette defined (1 "black", 2 "#d62728", 3 "#1f77b4", 4 "#2ca02c")
    set cbrange [0.5:4.5]
    set cbtics 1
    set cblabel "distinct atoms"
    p3 = "'trace_fc.dat' index 1 u 5:6:8 w p pt 6 ps 0.9 palette z title 'FC3'"
    p4 = "'trace_fc.dat' index 2 u 5:6:8 w p pt 8 ps 0.9 palette z title 'FC4'"
    if (n3 > 0 && n4 > 0) { eval "plot ".p3.", ".p4 } else { if (n3 > 0) { eval "plot ".p3 } else { eval "plot ".p4 } }
    pause -1
}
if (n4 > 0) {
    unset logscale y
    set format y "%g"
    set yzeroaxis
    set xlabel "Largest pair distance in the quadruplet (Ang)"
    set ylabel "sum_{ab} (Chi_{aabb}+Chi_{abab}+Chi_{abba})/3  (eV/Ang^4)"
    set title "FC4 isotropic part (5 x the average of Chi(n,n,n,n) over directions n)"
    plot 'trace_fc.dat' index 2 u 5:7:8 w p pt 8 ps 0.9 palette z notitle, 0 w l lt 5 lw 2 notitle
    pause -1
}
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
splot file0 u 2:3:4 with p pt 7 ps 0.2 lc rgb "yellow" title 'Atompos'   ,\
      file2 u 2:3:4 with p pt 7 ps 0.4 lc rgb "blue"   title 'R-vectors' ,\
      file3 u 2:3:4 with p pt 4 ps 0.3 lc rgb "green"  title 'SC atoms ' ,\
      file4 u 2:3:4 with p pt 2 ps 0.1 lc rgb "gray"   title 'WScell of SC'
pause -1
# one plot per primitive atom i: its FC2 bonds (springs<i>.dat: r_i - r_j) drawn from atom i to each
# neighbour j, and those neighbours. Atom i is put at the centre of the Wigner-Seitz cell of the supercell
# (the origin of WSR_boundary.xyz), not at its absolute position, so the neighbours r_j - r_i show the
# FC2 range against that cell (2026-11-03, user)
nat = int(system("ls springs[0-9]*.dat 2>/dev/null | wc -l"))
do for [i=1:nat] {
    fs = sprintf('springs%d.dat', i)
    xi = 0.0; yi = 0.0; zi = 0.0
    set title sprintf("FC2 neighbours of atom %d", i)
    splot fs u (xi):(yi):(zi):(-$4):(-$5):(-$6) with vectors nohead lt 1 title sprintf("%d-j", i) ,\
          '+' u (xi):(yi):(zi) every ::0::0 with p pt 7 ps 0.5 lc rgb "red" title sprintf('atom %d', i) ,\
          file2 u 2:3:4 with p pt 7 ps 0.3 lc rgb "blue" title 'R-vectors' ,\
          file4 u 2:3:4 with p pt 2 ps 0.1 lc rgb "gray" title 'WScell of SC'
#         fs u (xi-$4):(yi-$5):(zi-$6) with p pt 7 ps 0.3 lc rgb "black" title 'neighbours' ,\
#   pause -1
}
