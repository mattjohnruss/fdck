reset
set clip two
unset key

fn(i) = sprintf("output_%05i.csv", i)
fn_exact(i) = sprintf("output_exact_%05i.csv", i)

if(exists("t_inc")) {
}
else {
    t_inc = 1
}

splot for [i = 0:n:t_inc] fn(i) u 1:(i):j+1 w l lc j lw 1 ti columnheader
