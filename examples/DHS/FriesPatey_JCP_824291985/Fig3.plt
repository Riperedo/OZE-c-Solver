reset

set xlabel "r/{/Symbol s}"
set ylabel "h^{110}(r)"
set key top right

set xrange [1:4]
set yrange [-0.5:3.5]

plot "fig3_rho_0.8_mu2_2.0.dat" using 1:2 pt 7 ps 2 lc rgb "black" t "Monte Carlo", \
    "hdr_rho_0.8_mu2_2.0_LHNC.dat" using 1:3 w l lw 2 lc rgb "red" t "LHNC", \
    "hdr_rho_0.8_mu2_2.0_MSA.dat" using 1:3 w l lw 2 lc rgb "blue" t "MSA",