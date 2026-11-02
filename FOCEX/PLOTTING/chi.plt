set term pdfcairo enhanced font "Arial,16"
set output 'chi.pdf'
set key left

set xlabel "Frequency (1/cm)"
set grid xtics
set ylabel "Im({/Symbol c})"
 plot 'chi_imag.dat' u 1:2  w l ls 1 title '{/Symbol e}_{XX}' ,\
      'chi_imag.dat' u 1:6  w l ls 5 title '{/Symbol e}_{YY}' ,\
      'chi_imag.dat' u 1:10 w l ls 6 title '{/Symbol e}_{ZZ}' ,\
      'chi_imag.dat' u 1:11 w l ls 3 title 'n' ,\
      'chi_imag.dat' u 1:12 w l ls 4 title 'k' ,\
      'chi_imag.dat' u 1:3 w l ls 9 title 'OD' ,\
      'chi_imag.dat' u 1:4 w l ls 9 notitle,\
      'chi_imag.dat' u 1:5 w l ls 9 notitle,\
      'chi_imag.dat' u 1:7 w l ls 9 notitle,\
      'chi_imag.dat' u 1:8 w l ls 9 notitle,\
      'chi_imag.dat' u 1:9 w l ls 9 notitle
#pause -1

set ylabel "Re({/Symbol c})"
 plot 'chi_real.dat' u 1:2  w l ls 1 title '{/Symbol e}_{XX}' ,\
      'chi_real.dat' u 1:6  w l ls 5 title '{/Symbol e}_{YY}' ,\
      'chi_real.dat' u 1:10 w l ls 6 title '{/Symbol e}_{ZZ}' ,\
      'chi_real.dat' u 1:11 w l ls 3 title '{/Symbol e}_{normal}' ,\
      'chi_real.dat' u 1:12 w l ls 4 title 'R' ,\
      'chi_real.dat' u 1:3 w l ls 9 title 'OD' ,\
      'chi_real.dat' u 1:4 w l ls 9 notitle,\
      'chi_real.dat' u 1:5 w l ls 9 notitle,\
      'chi_real.dat' u 1:7 w l ls 9 notitle,\
      'chi_real.dat' u 1:8 w l ls 9 notitle,\
      'chi_real.dat' u 1:9 w l ls 9 notitle
#pause -1
