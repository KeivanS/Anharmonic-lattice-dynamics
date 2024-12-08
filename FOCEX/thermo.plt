set zeroaxis
set xlabel "strain(%)"
set ylabel "Free energy (kJ/mol)"
plot "temp_free.dat" u (100*$2):3 w p title "Free Energy",\
     "temp_free.dat" u (100*$2):4 w p title "Pressure"
replot "thermo_QHA.dat" u (100*$5):3 w lp notitle
pause -1
set xlabel "Temperature(K)"
set ylabel "Thermodynamic properties"
plot 'thermo_QHA.dat' u 1:2 title "E(kJ/mol)" ,\
     'thermo_QHA.dat' u 1:3 title "F(kJ/mol)"
pause -1
plot 'thermo_QHA.dat' u 1:6 title "Cv(J/K.mol)" ,\
     'thermo_QHA.dat' u 1:($7*1000) title "S(J/K.mol)"
pause -1
plot 'thermo_QHA.dat' u 1:4 title "Volume (Ang^3)"
pause -1
plot 'thermo_QHA.dat' u 1:($8*1E5) title "10^5 CVE(1/K)"
pause -1
plot 'thermo_QHA.dat' u 1:($9) title "Gruneisen"
pause -1
plot 'thermo_QHA.dat' u 1:($10) title "Bulk modulus(GPa)"
pause -1
