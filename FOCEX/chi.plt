set xlabel "Frequency (1/cm)"
set grid xtics
set ylabel "Imag(Xi)"
 plot 'chi_imag.dat' u 1:2 w l ls 1 title 'Diagonal' ,\
      'chi_imag.dat' u 1:6  w l ls 1 notitle ,\
      'chi_imag.dat' u 1:10 w l ls 1 notitle ,\
      'chi_imag.dat' u 1:3 w l ls 2 title 'off-diagonal' ,\
      'chi_imag.dat' u 1:4 w l ls 2 notitle,\
      'chi_imag.dat' u 1:5 w l ls 2 notitle,\
      'chi_imag.dat' u 1:7 w l ls 2 notitle,\
      'chi_imag.dat' u 1:8 w l ls 2 notitle,\
      'chi_imag.dat' u 1:9 w l ls 2 notitle
pause -1

set ylabel "Real(Xi)"
 plot 'chi_real.dat' u 1:2 w l ls 1 title 'Diagonal' ,\
      'chi_real.dat' u 1:6  w l ls 1 notitle ,\
      'chi_real.dat' u 1:10 w l ls 1 notitle ,\
      'chi_real.dat' u 1:3 w l ls 2 title 'off-diagonal' ,\
      'chi_real.dat' u 1:4 w l ls 2 notitle,\
      'chi_real.dat' u 1:5 w l ls 2 notitle,\
      'chi_real.dat' u 1:7 w l ls 2 notitle,\
      'chi_real.dat' u 1:8 w l ls 2 notitle,\
      'chi_real.dat' u 1:9 w l ls 2 notitle
pause -1
