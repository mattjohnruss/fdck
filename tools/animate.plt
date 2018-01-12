set term gif size 1280,720 animate delay 2

set output "animate.gif"
#unset key

set key top right

set yrange [0:2]

frames = system("ls -1 output*.csv | wc -l")
steady_frame = system("ls -1 output_steady*.csv | wc -l")
frames = frames - steady_frame

do for [i=0:frames-1] {
    set label 1 "".i."" at graph 1,1
    filename = sprintf("output_%05i.csv", i)
    plot filename u 1:2 w l ti columnhead, "" u 1:3 w l ti columnhead, "" u 1:4 w l ti columnhead, "" u 1:5 w l ti columnhead
}

set output
#!gwenview animate.gif
