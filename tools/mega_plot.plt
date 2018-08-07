set term gif size 1920, 1080 animate delay 3
set output "animate.gif"

set clip two

set key top right

prefixes = "const_chemotaxis_const_chemokinesis/ const_chemotaxis_hill_chemokinesis/ hill_chemotaxis_const_chemokinesis/ hill_chemotaxis_hill_chemokinesis/"

fn(prefix, i) = sprintf("%soutput_%05i.csv", prefix, i)
fn_exact(prefix, i) = sprintf("%soutput_exact_%05i.csv", prefix, i)
fn_steady(prefix) = sprintf("%soutput_steady.csv", prefix)

n_var = words(vars)

if(exists("n_inc")) {
}
else {
    n_inc = 1
}

do for [i = 0:n:n_inc] {
    set multiplot layout 2, 2

    set label 1 sprintf("%i",i) at graph 0.05,0.95

    do for [p = 1:words(prefixes)] {
        prefix = word(prefixes, p);
        
        title = system("echo ".prefix." | sed 's/_/ /g' | tr -d '\\\\/'")
        set title title

        #if(exist("steady")) {
            #plot for [j = 1:n_var] fn(prefix, i)     u 2:2+word(vars,j) w l lc j lw 2 ti columnheader, \
                 #for [j = 1:n_var] fn_steady(prefix) u 2:2+word(vars,j) w l lc j lw 2 dt 3 noti
        #}
        #else {
            #if(exist("exact")) {
                #plot for [j = 1:n_var] fn(prefix, i) u 2:2+word(vars,j) w l lc j lw 2 ti columnheader, \
                     #for [j = 1:n_var] fn(prefix, i) u 2:(exactf($2))   w l lc j dt 2 ti "exact"
            #}
            #else {
                #if(exist("exact_file")) {
                    #plot for [j = 1:n_var] fn(prefix, i)       u 2:2+word(vars,j) w l lc j lw 2 ti columnheader, \
                         #for [j = 1:n_var] fn_exact(prefix, i) u 2:2+word(vars,j) w l lc j lw 2 dt 2 ti columnheader
                #}
                #else {
                    plot for [j = 1:n_var] fn(prefix, i) u 2:2+word(vars,j) w l lc j lw 2 ti columnheader
                #}
            #}
        }

    unset multiplot
}
