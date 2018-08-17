#reset
set clip two
unset key

set xlabel "x"
set ylabel "timestep"

if(exist("title")) {
    set title title
}

fn(i) = sprintf("output_%05i.csv", i)
fn_exact(i) = sprintf("output_exact_%05i.csv", i)

if(exists("t_inc")) {
}
else {
    t_inc = 1
}

splot for [i = 0:n:t_inc] fn(i) u 2:(i):j+2 w l lc j lw 1 ti columnheader
