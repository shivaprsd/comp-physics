datafile = 'test.dat'
frames = "0 480 560 720 880 960 1040 1200 1520"
set term png size 1000,1280 font "Sans Serif"
set style line 10 lw 3 lc black
set border ls 10
set margin 0,0,0,0
set yrange [-2:2]
#set tics scale 0 font "Avenir,16" enhanced
unset tics

#stats datafile
do for [j in frames] {
    set out 'tmp/img'.sprintf("%.4d", int(j)).'.png'
    set label at graph 0.04, graph 0.94 font ",42" sprintf("%d", int(j))
    plot datafile i int(j) u 1 w l ls 10 not, '' i int(j) u 2 w l ls 10 not
    unset label
}
