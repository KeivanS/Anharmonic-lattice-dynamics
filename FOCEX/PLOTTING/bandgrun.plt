# Dispersion and mode Gruneisen parameters along the band path, sharing the k axis, one colour per
# branch in both panels (branch order and colours as bandgen.plt). Placeholders in capitals are filled
# by allplots; by hand: set the numbers (the path labels are read from KTICS.BS).
# bs_grun.dat: nk, dk (path length), k(1:3), then nb frequencies (1/cm), then nb Gruneisen parameters
set term pdfcairo size 5in,6in enhanced font "Arial,14"
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

# nice tic step for a range r (about 5 tics); the label at the join of the panels is left out
nice(r) = (ex = 10**floor(log10(r/5.)), fr = r/5./ex, (fr < 1.5 ? 1 : fr < 3 ? 2 : fr < 7 ? 5 : 10)*ex)
# path labels from KTICS.BS (set xtics (...)): on top of the bands (x2) and below the Gruneisen panel (x)
kt = system("cat KTICS.BS")
ktx2 = "set x2tics " . kt[strstrt(kt, "("):strlen(kt)]
fmax = MAX_FREQUENCY_IN_INVERSE_CM
gmin = BSGRUN_MIN
gmax = BSGRUN_MAX
xl = 0.16; xr = 0.84; ytop = 0.91; ymid = 0.52; ybot = 0.08

set multiplot title "Phonons of NAME"
# ---- bands (top) ----
set lmargin at screen xl; set rmargin at screen xr
set tmargin at screen ytop; set bmargin at screen ymid
set xrange [0 : MAX_K_VALUE ]
set x2range [0 : MAX_K_VALUE ]
unset xtics
eval(ktx2)
set grid xtics ytics x2tics lt 1 lw 0.6 lc rgb "#303030"   # k-path lines dark enough to see
set ylabel "Frequency (1/cm)"
set yrange [0:fmax]
set ytics nice(fmax), nice(fmax) nomirror          # no label at 0 (the join)
set y2range [0:fmax*3/100.]
set y2tics nice(fmax*3/100.), nice(fmax*3/100.) nomirror
set y2label "Frequency (THz)"
plot for [b=1:nb] 'bs_grun.dat' u 2:(column(5+b)) w l ls (1+(b-1)%8) notitle

# ---- Gruneisen parameters (bottom), touching the bands ----
set tmargin at screen ymid; set bmargin at screen ybot
unset x2tics
unset y2tics
unset y2label
eval(kt)
set grid xtics ytics x2tics lt 1 lw 0.6 lc rgb "#303030"   # k-path lines dark enough to see
set ylabel "Gruneisen parameter"
set yrange [gmin:gmax]
set ytics gmin, nice(gmax-gmin), gmax - 0.5*nice(gmax-gmin) mirror   # no label at gmax (the join)
set xzeroaxis
plot for [b=1:nb] 'bs_grun.dat' u 2:(column(5+nb+b)) w l ls (1+(b-1)%8) notitle
unset multiplot
