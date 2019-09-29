set term qt size 800,600 font "/Library/Fonts/Palatino.ttc,12"
set tics font "Palatino,12"
set key font "Palatino,14" height 2 spacing 1.5
set grid

f(x) = a * x + b
g(x) = c * x + d
h(x) = e * x + f

fit f(x) 'woutput.dat' u (log(2.268 - $1)):(log($2)) via a, b
fit g(x) 'woutput.dat' u (log(2.27 - $1)):(log($3)) via c, d
fit [-3.7:] h(x) 'woutput.dat' u (log($1 - 2.268)):(log($3)) via e, f

set xrange [-5.5:]
set title font "Palatino,24"
plot 'woutput.dat' u (log(2.27 - $1)):(log($3)) pt 6 lc black t "T < T_C", \
     [-5:] g(x) t sprintf("m = %.3f ", c) at beg, \
     'woutput.dat' u (log($1 - 2.268)):(log($3)) pt 7 lc black t "T > T_C", \
     [-5:] h(x) t sprintf("m = %.3f ", e) at beg
