set term pdfcairo enhanced font "Arial,16"
set output 'chi.pdf'
set key left

set xlabel "Frequency (1/cm)"
set grid xtics
set ylabel "Im({/Symbol c})"
# peak labels: the IR-active TO modes of chi_modes.dat (FOCEX, born_flag 2,3,6,7,8; the poles of chi): irrep,
# IR (+R if also Raman active), frequency; dashed line at each (2026-11-03, user). Raman-only modes give no peak.
if (system("test -f chi_modes.dat && echo 1") eq "1") {
   nm = int(system("awk '!/^#/ && $3==\"yes\"' chi_modes.dat | wc -l"))
   do for [k=1:nm] {
      ln = system(sprintf("awk '!/^#/ && $3==\"yes\"' chi_modes.dat | sed -n %dp", k))
      w = real(word(ln,1))
      lab = sprintf("%s (IR%s) %.0f", word(ln,7), (word(ln,4) eq "yes" ? "+R" : ""), w)
      set arrow k from w, graph 0 to w, graph 1 nohead dt 2 lc rgb "gray"
      set label k lab at w, graph 0.97 rotate by 90 right offset -0.7,0 font ",11"
   }
}
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
unset arrow
unset label

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
