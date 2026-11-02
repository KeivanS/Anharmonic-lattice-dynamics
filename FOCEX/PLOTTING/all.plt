# all placeholders are in capital letters and need to be edited
set term pdfcairo #postscript solid color enhanced # pdfcairo
set output "all-NAME.pdf"
set  multiplot layout 1,3 margins 0.10 ,0.95,0.15,0.90 spacing 0.13,0.02  title "Phonons of NAME" #, 0.01 , 0.01

# freq versus gruneisen ===========================

#set origin 0,0
#set size 0.3 ,1
#set lmargin 0
#set rmargin 0
a=MAX_FREQUENCY_IN_INVERSE_CM
set yrange [0:a]
set ylabel 'Frequency (1/cm)'
#set label  'Gruneisen' at graph 0.2,0.8 font ',8'
set ytics nomirror
set grid xtics
set grid ytics
set xtics in 1
set xlabel "Gruneisen"
set xrange [GRUNEISEN_MIN:GRUNEISEN_MAX]

# ibz_grun.dat: 5 columns (index, dk, k), then the nb frequencies, then the nb Gruneisen parameters;
# nb (= 3 x atoms) is read from the file
nb = int(system("awk '!/^ *#/ && NF>5 {print (NF-5)/2; exit}' ibz_grun.dat"))
plot for [b=1:nb] 'ibz_grun.dat' u (column(5+nb+b)):(column(5+b)) w p pt 1 ps 0.1 notitle

# fbz_dos is from tetrahedon method  =================================
# DOES NOT NEED EDITING AND SHOULD WORK BY DEFAULT
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
set y2range [0:a*3/100.]
set y2label 'Frequency (THz)'
set y2tics nomirror
set format y ""
############# REPLACE THE FOLLOWING LINE BY THE CONTENTS OF "KTICS.BS"  #################
 set xtics ( "L" 0.000,"G" 0.500,"K" 1.112,"X" 1.316,"W" 1.605,"L" 2.013,"K" 2.367,"W" 2.571,"X" 2.860,"G" 3.435 )
####### ACCORDINGLY REPLACE THE XRANGE BY THE LARGEST (LAST) NUMBER IN THE "KTICS.BS"
set xrange [0 : MAX_K_VALUE ]
set grid xtics
set xlabel "Wavenumbers"
set palette maxcolors 30  # adjust the number of colors as needed
set palette defined       # model RGB
set cbtics
#set bmargin 5
set colorbox user origin 0.92,0.05  # Places box at 95% right, 5% up from bottom
set colorbox size 0.02,0.85     
# coloured by the band Gruneisen parameters when bs_grun.dat exists (FC3 fitted), else by group velocity
# bs_grun.dat: nk, dk (path length), k(1:3), then nb frequencies (1/cm), then nb Gruneisen parameters
grun = int(system("test -s bs_grun.dat && echo 1 || echo 0"))
if (grun) {
    nbg = int(system("awk '!/^ *#/ && NF>5 {print (NF-5)/2; exit}' bs_grun.dat"))
    set cblabel "Gruneisen parameter"
    set cbrange [BSGRUN_MIN:BSGRUN_MAX]
    plot for [b=1:nbg] 'bs_grun.dat' u 2:(column(5+b)):(column(5+nbg+b)) w l lw 2.5 lc palette notitle
} else {
    set cblabel "Group Velocity(km/s)"
    set cbrange [0: MAX_GROUP_VELOCITY ] # same range as palette definition
    plot 'bs_freq.dat' u 3:7:11 w p pt 10 ps 0.23 lc palette z notitle
}

pause -1

unset multiplot
