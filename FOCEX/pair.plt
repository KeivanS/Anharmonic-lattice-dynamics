set xrange [1:13]
 set key top left
 set title 'Pair distances in ...'
 set xlabel 'Distance (Ang)'
 plot    'pairs.dat' u 1:($2) w l ls 2 title 'Atom1-'
 replot  'pairs.dat' u 1:(-$4) w l ls 3 title 'Atom2-'
pause -1
# set term postscript enhanced color font "arial,24"
# set output 'pairs.eps'
# replot

 set title 'FC2 versus Pair distances in ...'
 plot    'trace_fc.dat' u 1:($2) w p ps 2 notitle
pause -1
