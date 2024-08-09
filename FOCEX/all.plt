set term pdfcairo #postscript solid color enhanced # pdfcairo
set output "all-RbBr.pdf"
set  multiplot layout 1,3 margins 0.10 ,0.95,0.15,0.90 spacing 0.13,0.02  title "Phonons of RbBr" #, 0.01 , 0.01

# freq versus gruneisen ===========================

#set origin 0,0
#set size 0.3 ,1
#set lmargin 0
#set rmargin 0
a=150
set yrange [0:a]
set ylabel 'Frequency (1/cm)'
#set label  'Gruneisen' at graph 0.2,0.8 font ',8'
set ytics nomirror
set grid xtics
set grid ytics
set xtics in 1
set xlabel "Gruneisen"

plot 'ibz_grun.dat'   u 12:6  w p pt 1 ps 0.1 notitle ,\
     'ibz_grun.dat'   u 13:7  w p pt 1 ps 0.1 notitle ,\
     'ibz_grun.dat'   u 14:8  w p pt 1 ps 0.1 notitle ,\
     'ibz_grun.dat'   u 15:9  w p pt 1 ps 0.1 notitle ,\
     'ibz_grun.dat'   u 16:10 w p pt 1 ps 0.1 notitle ,\
     'ibz_grun.dat'   u 17:11 w p pt 1 ps 0.1 notitle

# fbz_dos is from tetrahedon method  =================================

#set origin 0.35,0
set lmargin at screen 0.3
set size ratio 3.0 #  0.2,1
unset ylabel
unset xrange
unset xtics
#set xrange [0:0.3]
set format y ""
set xlabel "DOS"
set key bottom
plot 'fbz_dos.dat' u 4:1  w l title "Total" ,\
     'fbz_dos.dat' u ($7+$6+$5):1  w filledcurve above y=0 title  "LA" ,\
     'fbz_dos.dat' u ($6+$5):1  w filledcurve  above y=0 title  "TA2" ,\
     'fbz_dos.dat' u 5:1  w filledcurve above y=0 title  "TA1"

#pause -1

# band structure  =================================

#set origin 0.56 ,0
#set size 0.80 , 1
set size ratio 1.2 #  0.2,1
set lmargin at screen 0.45
unset ylabel
set y2range [0:a*3/100]
set y2label 'Frequency (THz)'
set y2tics nomirror
set format y ""
set xrange [0 : 2.6 ]
 set xtics ( "G" 0.000,"X" 0.577,"K" 0.781,"G" 1.394,"L" 1.894,"W" 2.302,"X" 2.589 )
set grid xtics
set xlabel "Wavenumbers"
set palette maxcolors 30  # adjust the number of colors as needed
set palette defined       # model RGB
set cblabel "Group velovity(km/s)"
set cbtics
#set bmargin 5
set colorbox user origin 0.92,0.05  # Places box at 95% right, 5% up from bottom
set colorbox size 0.02,0.85     
set cbrange [0:6 ] # same range as palette definition

plot 'bands.dat' u 2:6:12  w l lc palette z  notitle ,\
     'bands.dat' u 2:7:13  w l lc palette z  notitle ,\
     'bands.dat' u 2:8:14  w l lc palette z  notitle ,\
     'bands.dat' u 2:9:15  w l lc palette z  notitle ,\
     'bands.dat' u 2:10:16 w l lc palette z  notitle ,\
     'bands.dat' u 2:11:17 w l lc palette z  notitle

pause -1

unset multiplot
