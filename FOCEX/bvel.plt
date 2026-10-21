set palette maxcolors 50  # adjust the number of colors as needed
set palette defined # model RGB

############# REPLACE THE FOLLOWING LINE BY THE CONTENTS OF "KTICS.BS"  #################
 set xtics ( "L" 0.000,"G" 0.500,"K" 1.112,"X" 1.316,"W" 1.605,"L" 2.013,"K" 2.367,"W" 2.571,"X" 2.860,"G" 3.435 )
####### ACCORDINGLY REPLACE THE XRANGE BY THE LARGEST (LAST) NUMBER IN THE "KTICS.BS"
set xrange [0 : MAX_K_VALUE ]


set grid xtics
set grid ytics

a=MAX_FREQUENCY_IN_INVERSE_CM
#set yrange [0:a]
set ylabel "Frequency (1/cm)"
set ytics nomirror
set y2range [0:a*3/100]
set y2label "Frequency (THz)"
set y2tics nomirror

set title "Phonon dispersion of NAME"

# to setup the colorbar
set cblabel "Group velovity(km/s)"
set cbtics
set cbrange [0:MAX_GROUP_VELOCITY] # same range as palette definition

plot 'bs_freq.dat' u 3:7:11 w p pt 7 ps 0.5 lc palette z notitle 
pause -1

set term pdfcairo font "Arial,16"
set output 'NAME.pdf'
replot

# ---------------- vgroup versus frequency ------------------
reset
set xlabel "Frequency (1/cm)"
set ylabel "Group Velovity (km/s)"
# THIS HAS AS MANY LINES AS THERE ARE BANDS; ACCORDINGLY COLUMN NYMBERS MUST BE EDITED
# First column goes from 6 to 6+3N-1 and second column goes from 6+3N to 6+6N-1 (N=#of atoms)
plot 'ibz_bands.dat' u ($7):($11/1) w p pt 7 ps 0.5 notitle ,\
       'bands.dat' u (sqrt($6)*521.11):12 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($7)*521.11):13 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($8)*521.11):14 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($9)*521.11):15 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($10)*521.11):16 w p pt 7 ps 0.3 notitle ,\
       'bands.dat' u (sqrt($11)*521.11):17 w p pt 7 ps 0.3 notitle


# pt=1(plus),2(cross),3(cross-plus),4(empty square),5 (diamond),6(triangle),7(circle),10(square)


# ---------------- vrgoup versus frequency ------------------

# pt=1(plus),2(cross),3(cross-plus),4(empty square),5 (diamond),6(triangle),7(circle),10(square)
