load "std.gp"
set size square
set fit quiet
set key bottom Left reverse height 5 autotitle columnhead
set angle degrees

set logscale y
set xrange [:16]
set xtics 2
set mxtics 2
set ytics format "10^{%L}"
set xlabel "Prediction time, T\n"
set ylabel "Normalized error, E\n\n"
set title "Time series: E vs T"

f1(x) = a1 * x + b1
f2(x) = a2 * x + b2
f3(x) = a3 * x + b3
f4(x) = a4 * x + b4
fit [:11] f1(x) ARG1 i 0 u ($0 + 1):(log($1)) via a1, b1
fit [:12] f2(x) ARG1 i 2 u ($0 + 1):(log($1)) via a2, b2
fit [:11] f3(x) ARG1 i 1 u ($0 + 1):(log($1)) via a3, b3
fit [:11] f4(x) ARG1 i 3 u ($0 + 1):(log($1)) via a4, b4

k = 16 / (9 * log(10))
set label at 1, 2e-4 rot by atan(k * a1) sprintf("m = %.3f", a1) font font.",14"
set label at 1, 5e-8 rot by atan(k * a2) sprintf("m = %.3f", a2) font font.",14"
set label at 1, 2e-3 rot by atan(k * a3) sprintf("m = %.3f", a3) font font.",14"
set label at 1, 5e-6 rot by atan(k * a4) sprintf("m = %.3f", a4) font font.",14"

plot ARG1 i 0 u ($0 + 1):1 pt 4 lc 7 lw 1.2, [:13] exp(f1(x)) lc 3 lw 1.2 not, \
     "" i 2 u ($0 + 1):1 pt 12 ps 1.2 lc 1 lw 1.2, [:13] exp(f2(x)) lc 2 lw 1.2 not, \
     "" i 1 u ($0 + 1):1 pt 6 lc 6 lw 1.2, [:13] exp(f3(x)) lc 8 lw 1.2 not, \
     "" i 3 u ($0 + 1):1 pt 8 ps 1.2 lc 8 lw 1.2, [:12] exp(f4(x)) lw 1.2 lc 4 not
