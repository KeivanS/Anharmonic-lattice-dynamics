bin(x,w)=w*floor(x/w)
# for a gaussian
bw=0.2   # bin width
plot 'gauss_0_1.dat' u (bin($1,bw)):(1.0) smooth freq with boxes
pause -1
plot 'dispvel0.dat'  u (bin($2+$3+$4,0.04)):(1.0) smooth freq with boxes
pause -1
plot 'dispvel0.dat'  u (bin($5+$6+$7,0.1)):(1.0) smooth freq with boxes
pause -1

set zeroaxis
#set xrange[1:1000]
  plot 'ene_lj.dat' u 2:3   title "KE"
replot 'ene_lj.dat' u 2:4   title "PE"
replot 'ene_lj.dat' u 2:5   title "TE"
replot 'ene_lj.dat' u 2:($6*1.5)   title "3kT/2"
replot 300*8.625e-5*1.5 
pause -1

plot 'ene_lj.dat' u 2:($10*1)   title "MSD"
replot sqrt(x)/4.216
#replot 2:(sqrt($2-1000)/4.216)
pause -1  

  plot 'cur_lj.dat' u 2:3   title "jx"
replot 'cur_lj.dat' u 2:4   title "jy"
replot 'cur_lj.dat' u 2:5   title "jz"
replot 'ene_lj.dat' u 2:7   title "x"
replot 'ene_lj.dat' u 2:8   title "y"
replot 'ene_lj.dat' u 2:9   title "z"
pause -1