# kappa.plt -- THERMACOND plots (2026-11-09), filled and run by allplots in a THERMACOND folder (NAME = material):
#   kappa_T-NAME.pdf               the thermal conductivity versus T: RTA, full solution, Wigner (P + C)
#   kappa_wigner_modes-NAME.pdf    Wigner kappa from the lowest bands, stacked: TA1, TA1+TA2, TA1+TA2+LA, all
#   kappa_cumul_mfp_T<T>-NAME.pdf  cumulative kappa versus mean free path, one file per temperature
#   kappa_cumul_freq_T<T>-NAME.pdf cumulative kappa versus frequency, one file per temperature
# Every curve is the average (xx+yy+zz)/3. Data (THERMACOND 2.2): kappa_summary.dat (T, RTA avg xx yy zz,
# full avg xx yy zz, Wigner avg xx yy zz), kappa_wigner_modes.dat (T, then avg xx yy zz of TA1, TA1+TA2,
# TA1+TA2+LA, all), kappa_cumul_mfp_T*.dat (mfp RTA, RTA avg xx yy zz, mfp eff, full avg xx yy zz, mfp RTA,
# Wigner avg xx yy zz), kappa_cumul_freq_T*.dat (omega 1/cm, THz, RTA avg xx yy zz, full ..., Wigner ...).
set term pdfcairo font 'Arial,12' size 5.5in,4in
set grid
set key top right

set output 'kappa_T-NAME.pdf'
set title 'NAME: thermal conductivity, (xx+yy+zz)/3'
set xlabel 'T (K)'
set ylabel 'kappa (W/mK)'
plot 'kappa_summary.dat' u 1:2 w lp pt 7 ps 0.6 lw 2 t 'RTA', \
     ''                  u 1:6 w lp pt 5 ps 0.6 lw 2 t 'full solution', \
     ''                  u 1:10 w lp pt 9 ps 0.6 lw 2 t 'Wigner (P + C)'

set output 'kappa_wigner_modes-NAME.pdf'
set title 'NAME: Wigner kappa summed over the lowest bands'
set style fill solid 0.85 noborder
plot 'kappa_wigner_modes.dat' u 1:14 w filledcurves x1 lc rgb '#c6dbef' t 'all bands (total)', \
     ''                       u 1:10 w filledcurves x1 lc rgb '#6baed6' t 'TA1 + TA2 + LA', \
     ''                       u 1:6  w filledcurves x1 lc rgb '#2171b5' t 'TA1 + TA2', \
     ''                       u 1:2  w filledcurves x1 lc rgb '#08306b' t 'TA1', \
     ''                       u 1:14 w l lw 1.5 lc rgb 'black' notitle
set style fill empty

# one PDF per temperature; the temperature is the part of the file name between '_T' and '.dat'
tlabel(f) = substr(f, strstrt(f, '_T') + 2, strlen(f) - 4)
set logscale x
set format x '10^{%T}'
set xlabel 'mean free path (nanometres)'
fm = system('ls kappa_cumul_mfp_T*.dat 2>/dev/null')
do for [f in fm] {
   t = tlabel(f)
   set output 'kappa_cumul_mfp_T'.t.'-NAME.pdf'
   set title 'NAME: cumulative kappa versus mean free path, T = '.t.' K'
   set key top left
   # from 0.001 nm: modes of (almost) zero velocity have mfp ~ 0 and carry nothing
   plot f u ($1 > 1e-3 ? $1 : 1/0):2   w l lw 2 t 'RTA', \
        f u ($6 > 1e-3 ? $6 : 1/0):7   w l lw 2 t 'full solution (tau_{eff})', \
        f u ($11 > 1e-3 ? $11 : 1/0):12 w l lw 2 t 'Wigner (P + C)'
}
unset logscale x
set format x '% h'
set xlabel 'frequency (1/cm)'
ff = system('ls kappa_cumul_freq_T*.dat 2>/dev/null')
do for [f in ff] {
   t = tlabel(f)
   set output 'kappa_cumul_freq_T'.t.'-NAME.pdf'
   set title 'NAME: cumulative kappa versus frequency, T = '.t.' K'
   set key bottom right
   plot f u 1:3  w l lw 2 t 'RTA', \
        f u 1:7  w l lw 2 t 'full solution', \
        f u 1:11 w l lw 2 t 'Wigner (P + C)'
}
unset output
