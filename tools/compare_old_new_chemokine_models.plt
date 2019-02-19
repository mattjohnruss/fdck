file(prefix, i) = sprintf("%s/output_%05i", prefix, i)

# c_u  c_b  c_s  phi_i  phi_m  phi_c_u  phi_c_b
# 3    4    5    6      7
# 8    9    10   11     12     13       14

plot "< exec bash -c \"paste -d' ' res_ai/output_00100.csv <(cut -d' ' -f3- res_mb/output_00100.csv)\"" u 2:($8-$3) w p lw 2 ti col, \
     "" u 2:($9-$4) w p lw 2 ti col#, \
     #"" u 2:($10-$5) w p lw 2 ti col, \
     #"" u 2:($11-$6) w p lw 2 ti col, \
     #"" u 2:($12-$7) w p lw 2 ti col

#command(i) = sprintf("< exec bash -c \\\"paste -d' ' res_ai/output_%05i.csv <(cut -d' ' -f3- res_mb/output_%05i.csv)\\\"", i)

#plot command(100) u 2:($8-$3) w l lw 2 ti col, \
     #"" u 2:($9-$4) w l lw 2 ti col, \
     #"" u 2:($10-$5) w l lw 2 ti col, \
     #"" u 2:($11-$6) w l lw 2 ti col, \
     #"" u 2:($12-$7) w l lw 2 ti col
