if (!exists("fin"))    fin    = 'bands.dat'   # default is bands.dat
# number of bands: given with -e "nbands=N", else read from the file (5 columns, then nbands omega^2, nbands |v|)
if (!exists("nbands")) nbands = int(system(sprintf("awk '!/^ *#/ && NF>5 {print (NF-5)/2; exit}' %s", fin)))
if (!exists("fout"))   fout   = 'bands.png'   # default output is bands.png
print "nbands =", nbands
print "fin  =", fin
print "fout =", fout

col_start = 6
col_end = col_start + nbands - 1
print "range  =", col_start, col_end

set key outside
set grid
set xlabel "Q"
set ylabel "Eigenvalues"


# --- styles (ADD THIS BLOCK) ---
set style data lines
set style line 1 lc rgb "#E41A1C" lt 1 lw 3     # red
set style line 2 lc rgb "#377EB8" lt 1 lw 3     # blue
set style line 3 lc rgb "#4DAF4A" lt 1 lw 3     # green
set style line 4 lc rgb "#984EA3" lt 1 lw 3     # purple
set style line 5 lc rgb "#FF7F00" lt 1 lw 3     # orange
set style line 6 lc rgb "#999999" lt 1 lw 3     # gray
set style line 7 lc rgb "#F781BF" lt 1 lw 3     # pink
set style line 8 lc rgb "#A65628" lt 1 lw 3     # brown
# -------------------------------

### Build the plot command
lc = 0
plotcmd = ""
p2      = ""

do for [c = col_start:col_end] {
    lc = lc + 1

    # wrap around if more bands than defined styles
    ls_idx = 1 + (lc-1) % 8

    entry  = sprintf("'%s' using 2:%d with l ls %d notitle",          fin, c, ls_idx)
    entry2 = sprintf("'%s' using 2:(521.1*sgn($%d)*sqrt(abs($%d))) with l ls %d notitle", fin, c, c, ls_idx)

    if (c == col_start) {
        plotcmd = plotcmd . entry
        p2      = p2      . entry2
    } else {
        plotcmd = plotcmd . "," . entry
        p2      = p2      . "," . entry2
    }
}

# optional bounds: kmax (last value of KTICS.BS), fmax (1/cm; the omega^2 axis gets (fmax/521.1)^2)
if (exists("kmax")) set xrange [0:kmax]
if (exists("fmax")) set yrange [*:(fmax/521.1)**2]
eval("plot " . plotcmd)
pause -1

### Set output and style (fmt='pdf' for a pdf, as allplots does; default png)
if (!exists("fmt")) fmt = 'png'
if (fmt eq 'pdf') {
    set terminal pdfcairo size 6in,4in enhanced font ",11"
} else {
    set terminal pngcairo size 1200,800 enhanced font ",16"
}
set output sprintf("%s", fout)
replot

set ylabel "Frequency (1/cm)"
if (exists("fmax")) set yrange [*:fmax]
set y2tics
set link y2 via y/100.0*3.0 inverse y*100/3.
set y2label "Frequency (THz)"
set output sprintf("freq_%s", fout)
eval("plot " . p2 . " axes x1y1")

