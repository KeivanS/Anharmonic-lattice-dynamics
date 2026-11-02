set term pdfcairo
set output 'qha.pdf'

set multiplot layout 1,3 margins 0.05 ,0.98 ,0.15,0.90 spacing 0.05,0.02  

set xlabel 'Temperature(K)'
set xtics 250

#set size ratio 2
#set ylabel 'CTE*10^6'
#set yrange [:3]
plot 'thermo_QHA.dat' u 1:(100000*$9) lw 2 title "10^5 CTE" ,\
     'thermo_QHA.dat' u 1:10 lw 2 title "Gruneisen " ,\
     'thermo_QHA.dat' u 1:($11) lw 2 title "B_T(T)/B_T(0K)"
#pause -1

#set size ratio 2
unset yrange
#set ylabel 'Strain'
plot 'thermo_QHA.dat' u 1:(100*($4-1)) lw 2 title "dV/V(%)" ,\
     'thermo_QHA.dat' u 1:(300*$5)  w p ps 0.5 title "3*Strain(%)"
#pause -1

#set size ratio 2
#set ylabel 'Energies'
plot 'thermo_QHA.dat' u 1:2 lw 2 title "E(kJ/mol)" ,\
     'thermo_QHA.dat' u 1:3 lw 2 title "F(kJ/mol)" ,\
     'thermo_QHA.dat' u 1:6 lw 2 title "C_V(J/mol.K)" ,\
     'thermo_QHA.dat' u 1:7 lw 2 title "C_P(J/mol.K)" ,\
     'thermo_QHA.dat' u 1:($8) lw 2 title "S(J/mol.K)"
#pause -1

unset multiplot

reset

set yrange [-30:20]
set xlabel 'Strain (%)'
set ylabel 'Free Energy (KJ/mol)'
set key top center
# temp_free.dat has one block per temperature without blank lines between them; this filter adds them
# (gnuplot needs them for the contour map, and to draw each F(V), P(V) curve separately)
tf = '< awk ''{ if ($1 !~ /^#/ && NF > 0) { if (seen && $1 != p) print ""; p = $1; seen = 1 } print }'' temp_free.dat'
plot tf u (100*$2):3 w l title 'F(V)' ,\
     tf u (100*$2):4 w l title 'P(V)'  ,\
     'thermo_QHA.dat' u (100*$5):3 w l lw 2 title "F_{eq}(kJ/mol)" ,\
  0 w l lt 8 lw 2 notitle
#pause -1

unset style data
set xlabel "Temperature (K)"
set ylabel "Strain (\%)"
set title "Free energy contours"
unset yrange
set pm3d map
set palette defined
set view map
set contour surface   
set cntrparam levels auto 10
unset clabel                       # all contour lines in one style
set autoscale xfix
set autoscale yfix
set key top right textcolor rgb "white"
splot tf u 1:(100*$2):3 w pm3d notitle ,\
      tf u 1:(100*$2):3 w l nosurface lc rgb "black" notitle ,\
      'thermo_QHA.dat' u 1:(100*$5):3 w l nocontours lw 2 lc rgb "white" title "equilibrium strain"
#pause -1

