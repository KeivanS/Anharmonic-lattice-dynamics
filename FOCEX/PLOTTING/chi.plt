set term pdfcairo enhanced font "Arial,16"
set output 'chi.pdf'
set key left

set xlabel "Frequency (1/cm)"
set grid xtics
set ylabel "Im({/Symbol c})"
# peak labels: the IR-active TO modes of chi_modes.dat (FOCEX, born_flag 2,3,6,7,8; the poles of chi): irrep,
# IR (+R if also Raman active), frequency; dashed line at each (2026-11-03, user). Raman-only modes give no peak.
# Labels closer than dmin (about 1.5 label heights, 6 % of the frequency range of chi_imag.dat) to the
# previous one are moved to the right by dmin, so that they do not overlap; the dashed line stays at the
# frequency (2026-11-04, user). Modes are taken in increasing frequency.
if (system("test -f chi_modes.dat && echo 1") eq "1") {
   nm = int(system("awk '!/^#/ && $3==\"yes\"' chi_modes.dat | wc -l"))
   xlo = real(system("awk '$1!~/^#/ && NF>1 {x=$1+0; if(n==0||x<m)m=x; n++} END{print m}' chi_imag.dat"))
   xhi = real(system("awk '$1!~/^#/ && NF>1 {x=$1+0; if(n==0||x>m)m=x; n++} END{print m}' chi_imag.dat"))
   dmin = 0.06*(xhi-xlo)
   xl = -1e30
   do for [k=1:nm] {
      ln = system(sprintf("awk '!/^#/ && $3==\"yes\"' chi_modes.dat | sort -g | sed -n %dp", k))
      w = real(word(ln,1))
      lab = sprintf("%s (IR%s) %.0f", word(ln,7), (word(ln,4) eq "yes" ? "+R" : ""), w)
      xl = (w-xl < dmin) ? xl+dmin : w
      set arrow k from w, graph 0 to w, graph 1 nohead dt 2 lc rgb "gray"
      set label k lab at xl, graph 0.97 rotate by 90 right offset -0.7,0 font ",11"
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
