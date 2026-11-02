# Dispersion and mode Gruneisen parameters along the band path, sharing the k axis, one colour per
# branch in both panels (branch order and colours as bandgen.plt). Placeholders in capitals are filled
# by allplots; by hand: paste KTICS.BS over the xtics line and set the numbers.
# bs_grun.dat: nk, dk (path length), k(1:3), then nb frequencies (1/cm), then nb Gruneisen parameters
set term pdfcairo size 6in,6.5in enhanced font ",11"
set output 'bandgrun-NAME.pdf'
nb = int(system("awk '!/^ *#/ && NF>5 {print (NF-5)/2; exit}' bs_grun.dat"))

set style line 1 lc rgb "#E41A1C" lt 1 lw 2     # red
set style line 2 lc rgb "#377EB8" lt 1 lw 2     # blue
set style line 3 lc rgb "#4DAF4A" lt 1 lw 2     # green
set style line 4 lc rgb "#984EA3" lt 1 lw 2     # purple
set style line 5 lc rgb "#FF7F00" lt 1 lw 2     # orange
set style line 6 lc rgb "#999999" lt 1 lw 2     # gray
set style line 7 lc rgb "#F781BF" lt 1 lw 2     # pink
set style line 8 lc rgb "#A65628" lt 1 lw 2     # brown

set multiplot layout 2,1 margins 0.12,0.88,0.08,0.94 spacing 0,0.02 title "Phonons of NAME"
############# REPLACE THE FOLLOWING LINE BY THE CONTENTS OF "KTICS.BS"  #################
 set xtics ( "L" 0.000,"G" 0.500,"K" 1.112,"X" 1.316,"W" 1.605,"L" 2.013,"K" 2.367,"W" 2.571,"X" 2.860,"G" 3.435 )
set xrange [0 : MAX_K_VALUE ]
set grid xtics ytics
set format x ""
set ylabel "Frequency (1/cm)"
set yrange [0:MAX_FREQUENCY_IN_INVERSE_CM]
set y2range [0:MAX_FREQUENCY_IN_INVERSE_CM*3/100.]
set y2tics nomirror
set ytics nomirror
set y2label "Frequency (THz)"
plot for [b=1:nb] 'bs_grun.dat' u 2:(column(5+b)) w l ls (1+(b-1)%8) notitle

set format x
unset y2tics
unset y2label
set ytics mirror
set ylabel "Gruneisen parameter"
set yrange [BSGRUN_MIN:BSGRUN_MAX]
set xzeroaxis
plot for [b=1:nb] 'bs_grun.dat' u 2:(column(5+nb+b)) w l ls (1+(b-1)%8) notitle
unset multiplot
