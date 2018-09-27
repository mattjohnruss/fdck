reset

stats 'trace.dat' skip 1
#n = floor(sqrt(STATS_columns))

set multiplot layout 3, 3
#set multiplot layout 2, 2

do for [i=2:STATS_columns] {
    plot "trace.dat" u 1:i w l lw 1.5 ti col
}

unset multiplot
