set term qt size 800,600 font "/Library/Fonts/Palatino.ttc,12"
set tics font "Palatino,12"
set key font "Palatino,14" height 1.2 spacing 1.2 left
set fit quiet
set style line 1 lc black lw 1.5
set size square
set grid
set title "Ising 2D critical exponent: β" font "Palatino,18"
set xlabel "log(T_c - T)" font "Palatino,16"
set ylabel "log M" font "Palatino,16"
set label "T_c = 2.27" font "Palatino,14" at graph 1, graph 0 offset -1.2, 1.2 right

f1(x) = b1 * x + p1
f2(x) = b2 * x + p2
f3(x) = b3 * x + p3
f4(x) = b4 * x + p4
xmin = -4       #-1.7
xmax = -1.5     #-0.7

fit [xmin:xmax] f1(x) 'exp32.dat' i 1 via b1, p1
fit [xmin:xmax] f2(x) 'exp64.dat' i 1 via b2, p2
fit [xmin:xmax] f3(x) 'exp128.dat' i 1 via b3, p3
fit f4(x) 'exp256.dat' i 1 via b4, p4

plot [-5:] 'exp32.dat' i 1 pt 4 t "32×32", \
     'exp64.dat' i 1 pt 5 t "64×64", \
     'exp128.dat' i 1 lc black pt 6 t "128×128", \
     'exp256.dat' i 1 pt 7 t "256×256", \
     [-4:] f1(x) lw 1.5 t "β = –0.109  " at beg right, \
     [-4.4:-0.5] f2(x) lw 1.5 t "β = –0.118   " at end right, \
     [-4.6:-1] f3(x) lw 1.5 t "β = –0.117   " at end right, \
     [-4.8:-1.5] f4(x) lw 1.5 t "β = –0.118   " at end right
