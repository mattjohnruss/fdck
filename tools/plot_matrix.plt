reset

set xrange [-0.5:n*m-1+0.5]
set yrange [-0.5:n*m-1+0.5]

do for [i = 1:m-1] {
    pos = i*n - 0.5
    set arrow i from pos,-0.5 to pos,(n*m-1+0.5) nohead
    set arrow m-1+i from -0.5,pos to (n*m-1+0.5),pos nohead
}

unset key
set size square
#filename = "jac_t=0.100000_0.dat"
#load pal("parula")

if(exists("ps")) {
    plot filename w p pt 5 ps ps palette
}
else {
    plot filename w p pt 5 ps 1 palette
}
