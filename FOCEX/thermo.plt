set term pdfcairo font "Arial,16"
set output "thermo.pdf"

set zeroaxis
set xlabel "strain(%)"
set ylabel "Free energy (kJ/mol)"
set yrange [-20:50]
set key center top
plot "temp_free.dat" u (100*$2):3 w l lw 2 title "Free Energy",\
     "temp_free.dat" u (100*$2):4 w l lw 1 title "Pressure"
     "thermo_QHA.dat" u (100*$5):3 w p ps 0.5 pt 10 title "Equilibrium strain"
pause -1
set xlabel "Temperature(K)"
set ylabel "Thermodynamic properties"
unset yrange
set key top left
plot 'thermo_QHA.dat' u 1:2 lw 2 title "E(kJ/mol)" ,\
     'thermo_QHA.dat' u 1:3 lw 2 title "F(kJ/mol)" ,\
     'thermo_QHA.dat' u 1:6 lw 2 title "C_V(J/K.mol)" ,\
     'thermo_QHA.dat' u 1:($7) lw 2 title "C_P(J/K.mol)" ,\
     'thermo_QHA.dat' u 1:($8) lt 8 lw 2 title "S(J/K.mol)"
#pause -1
plot 'thermo_QHA.dat' u 1:($4-1)*100 lw 2 title "dV/V (%)" ,\
     'thermo_QHA.dat' u 1:($5-0)*300 w p ps 0.8 title "3 x strain (%)"
#pause -1
set key top right
plot 'thermo_QHA.dat' u 1:($9*1E5) lw 2 title "10^5 CVE(1/K)" ,\
     'thermo_QHA.dat' u 1:($10)  lw 2 title "Gruneisen" ,\
     'thermo_QHA.dat' u 1:($11) lw 2 title "B(T)/B_0"
#pause -1
