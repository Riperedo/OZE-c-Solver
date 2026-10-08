# plot_correlations.gp
# Comparison of g(r) and S(k) verifying exact numerical consistency

set terminal pdfcairo enhanced color font "Helvetica,10" size 6.0,4.5
set output "plots/fig_correlations_validation.pdf"

set multiplot layout 2,1 title "Numerical Consistency Verification: Cold vs Warm-Start ({/Symbol f} = 0.55)" font ",12"

# Top: Radial Distribution Function g(r)
set bmargin 2
set tmargin 2
set lmargin 8
set rmargin 3
set xlabel "r / {/Symbol s}"
set ylabel "g(r)"
set xrange [0.8:4.0]
set yrange [0.0:4.0]
set grid lc rgb "#e8e8e8" dt 2
set key top right

plot "data/gr_cold_phi_0.55.dat" using 1:2 with lines lw 3 lc rgb "#377eb8" title "Cold-Start Solution", \
     "data/gr_warm_phi_0.55.dat" using 1:2 with points pt 6 ps 0.7 lc rgb "#e41a1c" title "Warm-Start Solution"

# Bottom: Structure Factor S(k) vs Exact Analytical Wertheim Solution
set xlabel "k {/Symbol s}"
set ylabel "S(k)"
set xrange [0.0:25.0]
set yrange [0.0:4.2]
set key top right

plot "data/sk_analytical_phi_0.55.dat" using 1:2 with lines lw 2.5 lc rgb "#000000" title "Exact Wertheim-Thiele Analytical", \
     "data/sk_warm_phi_0.55.dat" using 1:2 with points pt 7 ps 0.6 lc rgb "#4daf4a" title "Warm-Start OZE Solver"

unset multiplot

set terminal pngcairo enhanced color font "Helvetica,10" size 1000,750
set output "plots/fig_correlations_validation.png"
replot
