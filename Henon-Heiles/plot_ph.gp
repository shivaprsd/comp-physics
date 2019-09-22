E = 1 / 12.0
i = 4

set term qt size 800,600 font "/Library/Fonts/Avenir.ttc,12"
set tics font "Avenir,12"
set key font "Avenir,12"
set size ratio -1
set grid

#set title "Phase space" font "Avenir,18"
set label "y" font "Avenir,16" at graph 1 offset char -1.5, -0.8
set label "ẏ" font "Avenir,16" at graph 0, graph 1 offset char -2.5, -0.5
set label sprintf("E = %.6f", E) font "Avenir,14" at graph 1, graph 1 offset char -14.5, -1.5
set xrange [-0.49:0.59]
set yrange [-0.49:0.49]
set xtics 0.1 scale 0.5
set ytics 0.1 scale 0.5
#set bmargin 4

plot sprintf("fig%d.dat", i) u 4:6 every 9 i 6 w lp pt 7 ps 0.2 lc black not

#unset label
#set label "x" font "Avenir,16" at graph 1 offset char -1.5, -0.8
#set label "y" font "Avenir,16" at graph 0, graph 1 offset char -2.5, -0.5
#set label "a = 1.6" font "Avenir,14" at graph 1, graph 1 offset char -8.5, -1.5
#set xrange [-0.99:0.99]
#set yrange [-0.99:0.99]
#set xtics 0.2
#set ytics 0.2
#set size ratio -1
#plot 'fig8.dat' pt 7 ps 0.3 lc black not

#set xtics 0.02
#set xrange [0:0.18]
#set yrange [0:1.05]
#set size noratio
#f(x) = a * x + b
#g(x) = (x >= 0.11 ? f(x) : 1)
#h(x) = (x <= 0.167 ? g(x) : 1/0)
#unset label
#set xlabel "ENERGY" font "Avenir,14"
#set ylabel "RELATIVE AREA" font "Avenir Bold,14" enhanced
#fit f(x) 'fig7.dat' u 1:($1 > 0.1 ? $2 : 1/0) via a,b
#plot 'fig7.dat' pt 7 lw 1.2 lc black not, h(x) lw 1.2 lc black not
