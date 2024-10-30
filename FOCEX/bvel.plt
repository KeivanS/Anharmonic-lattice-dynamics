set palette maxcolors 20  # adjust the number of colors as needed
set palette defined # model RGB

 set xrange [0:3.1]
 set xtics ( "G" 0.000,"K" 0.612,"X" 0.816,"G" 1.394,"L" 1.894,"X" 2.394,"W" 2.683,"L" 3.089 )
set grid xtics
set grid ytics

a=900
#set yrange [0:a]
set ylabel "Frequency (1/cm)"
set ytics nomirror
set y2range [0:a*3/100]
set y2label "Frequency (THz)"
set y2tics nomirror

set title "Phonon dispersion of MgO"

# to setup the colorbar
set cblabel "Group velovity(km/s)"
set cbtics
set cbrange [0:15] # same range as palette definition

plot 'bs_freq.dat' u 3:7:11 w p pt 7 ps 0.5 lc palette z notitle ,\
     'bs_freq.0'   u 3:7:11 w p pt 1 ps 0.2 lc palette z notitle 
pause -1

set term pdfcairo font "Arial,16"
set output 'MgO.pdf'
replot

# ---------------- vrgoup versus frequency ------------------
reset
set xlabel "Frequency (1/cm)"
set ylabel "Group Velovity (km/s)"
plot 'ibz_bands.dat' u ($7):($11/1) w p pt 7 ps 0.5 notitle ,\
       'bands.dat' u (sqrt($6)*521.11):12 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($7)*521.11):13 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($8)*521.11):14 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($9)*521.11):15 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($10)*521.11):16 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($11)*521.11):17 w p pt 7 ps 0.3 notitle

# 'veltest.dat' u (sqrt($2)*521.11):($6/1000) w p pt 7 ps 0.5 notitle ,\

# pt=1(plus),2(cross),3(cross-plus),4(empty square),5 (diamond),6(triangle),7(circle),10(square)


# ---------------- vrgoup versus frequency ------------------

# pt=1(plus),2(cross),3(cross-plus),4(empty square),5 (diamond),6(triangle),7(circle),10(square)
