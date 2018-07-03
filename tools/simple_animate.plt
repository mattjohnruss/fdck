set clip two

set key top right

fn(i) = sprintf("output_%05i.csv", i)
fn_exact(i) = sprintf("output_exact_%05i.csv", i)

n_var = words(vars)

do for [i = 0:n] {
    set label 1 sprintf("%i",i) at graph 0.05,0.95

    if(exist("steady")) {
        plot for [j = 1:n_var] fn(i)               u 1:1+word(vars,j) w l lc j lw 2 ti columnheader, \
             for [j = 1:n_var] "output_steady.csv" u 1:1+word(vars,j) w l lc j dt 2 noti
    }
    else {
        if(exist("exact")) {
            plot for [j = 1:n_var] fn(i)           u 1:1+word(vars,j) w l lc j lw 2 ti columnheader, \
                 for [j = 1:n_var] exactf(x)                          w l lc j dt 2 ti "exact"
        }
        else {
            if(exist("exact_file")) {
                plot for [j = 1:n_var] fn(i)           u 1:1+word(vars,j) w l lc j lw 2 ti columnheader, \
                     for [j = 1:n_var] fn_exact(i)     u 1:1+word(vars,j) w l lc j lw 2 dt 2 ti columnheader
            }
            else {
                plot for [j = 1:n_var] fn(i)           u 1:1+word(vars,j) w l lc j lw 2 ti columnheader
            }
        }
    }

    pause(0.02)
}
