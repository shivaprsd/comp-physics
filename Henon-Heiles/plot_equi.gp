set term qt size 800,600 font "/Library/Fonts/Avenir.ttc,12"
set tics font "Avenir,12"
set key font "Avenir,12"
set size ratio -1
set grid

set label "x" font "Avenir,16" at graph 1 offset char -1.5, -0.8
set label "y" font "Avenir,16" at graph 0, graph 1 offset char -2.5, -0.5
set label sprintf("%.4f", 1 / 100.0) font "Avenir,12" at first 0, first 0.18 center
set label sprintf("%.4f", 1 / 24.0) font "Avenir,12" at first 0, first 0.355 center
set label sprintf("%.4f", 1 / 12.0) font "Avenir,12" at first 0, first 0.535 center
set label sprintf("%.4f", 1 / 8.0) font "Avenir,12" at first 0, first 0.705 center
set label sprintf("%.4f", 1 / 6.0) font "Avenir,12" at first 0, first 1.04 center
set xrange [-0.99:0.99]
set yrange [-0.7:1.1]
set xtics 0.2 scale 0.5
set ytics 0.2 scale 0.5
set size ratio -1
plot 'fig2.dat' w l lw 1.2 lc black not, '' u (-$1):2 w l lw 1.2 lc black not
