splot "KPOINT.IBZ" u 2:3:4 w p ps 0.8 title "IBZ"  ,\
   "KPOINT_tet.MP" u 2:3:4 w p ps 0.4 title "FBZ" ,\
   "WSG_boundary.xyz" u 2:3:4 w p pt 2 ps 0.2 title "WS-cell"
pause -1
splot "KPOINT.IBZ" u 6:7:8 w p ps 0.8 title "IBZ"  ,\
   "KPOINT_tet.MP" u 6:7:8 w p ps 0.4 title "FBZ"
pause -1
