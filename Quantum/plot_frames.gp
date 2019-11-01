datafile = 'test.dat'
frames = "0 300 350 400 500 600 700 800 950"
set term png size 1000,1280 font "Sans Serif"
set style line 10 lw 4 lc black
set border ls 10
set margin 0,0,0,0
set yrange [-1:2]
set xrange [180:1234]
#set tics scale 0 font "Avenir,16" enhanced
unset tics

do for [j in frames] {
    set out 'tmp/img'.sprintf("%.4d", int(j)).'.png'
    #set label at graph 0.96, graph 0.94 right font ",48" "E = V_0"
    set label at graph 0.04, graph 0.94 font ",48" sprintf("%d", int(j))
    plot datafile i int(j + 1) w l ls 10 not, '' i 0 u ($1 / 1.5) w l ls 10 not
    unset label
}
