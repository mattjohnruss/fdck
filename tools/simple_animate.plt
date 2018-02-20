set clip two

fn(i) = sprintf("output_%05i.csv", i)

n_var = words(vars)

do for [i = 0:n] {
    set label 1 sprintf("%i",i) at graph 0.05,0.95

    if(exist("steady")) {
        plot for [j = 1:n_var] fn(i)               u 1:1+word(vars,j) w l ti columnheader, \
             for [j = 1:n_var] "output_steady.csv" u 1:1+word(vars,j) w l dt 2 noti
    }
    else {
        if(exist("exact")) {
            plot for [j = 1:n_var] fn(i)           u 1:1+word(vars,j) w l ti columnheader, \
                 for [j = 1:n_var] exactf(x)                          w l dt 2 ti "exact"
        }
        else {
            plot for [j = 1:n_var] fn(i)           u 1:1+word(vars,j) w l ti columnheader
        }
    }

    pause(0.02)
}
