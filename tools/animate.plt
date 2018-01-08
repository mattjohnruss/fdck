set term gif size 1280,720 animate delay 2

set output "animate.gif"
#unset key

set key top right

set yrange [0:1]

frames=system("ls -1 output*.dat | wc -l") - 1

do for [i=0:frames] {
    set label 1 "".i."" at graph 1,1
    plot "output_".i.".dat" u 1:2 w l ti "C_u", "" u 1:3 w l ti "C_b", "" u 1:4 w l ti "C_s", "" u 1:5 w l ti "φ"
}

set output
!gwenview animate.gif
