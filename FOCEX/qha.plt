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
     'thermo_QHA.dat' u 1:($11) lw 2 title "B(T)/B(0K)"
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
set key top center
plot 'temp_free.dat' u (100*$2):3 w l title 'F(V)' ,\
     'temp_free.dat' u (100*$2):4 w l title 'P(V)'  ,\
     'thermo_QHA.dat' u (100*$5):3 w l lw 2 title "F_{eq}(kJ/mol)" ,\
  0 w l lt 8 lw 2 notitle
#pause -1

unset style data
set xlabel "Temperature (K)"
set ylabel "Strain (\%)"
set title "Free energy contours"
set yrange [-0.05:0.10]
set pm3d map
set palette defined
set view map
set contour surface   
set cntrparam levels auto 10
splot 'temp_free.dat' u 1:2:3  w pm3d notitle
#pause -1

