load "parula_lines.pal"

fn(i) = sprintf("output_%05i.csv", i)

do for [i = n_min:n_max] {
    stats [*:*] [*:*] fn(i) u 1 skip 1 nooutput
    set label 1 sprintf("t = %.3f", STATS_min) at graph 0.05,0.95

    splot sprintf("output_%05u.csv", i) u 2:3:5 w surf ti col; pause 0.1
