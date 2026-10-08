# plot_gamma.gp
# Indirect correlation function gamma(r) across volume fractions phi in {0.50, 0.55, 0.60, 0.64}

set terminal pdfcairo enhanced color font "Helvetica,11" size 5.0,3.5
set output "plots/fig_gamma_profiles.pdf"

set title "Indirect Correlation Function {/Symbol g}(r) Profiles\n(Percus-Yevick Hard Spheres, N_{nodes}=4096)" font ",11"
set xlabel "Dimensionless Radial Distance r / {/Symbol s}"
set ylabel "Indirect Correlation Function {/Symbol g}(r)"
set grid lc rgb "#e8e8e8" dt 2

set xrange [0.0:4.5]
set yrange [-1.5:80.0]
set key top right box opaque

plot "data/gamma_phi_0.50.dat" using 1:2 with lines lw 2 lc rgb "#1b9e77" title "{/Symbol f} = 0.50", \
     "data/gamma_phi_0.55.dat" using 1:2 with lines lw 2 lc rgb "#d95f02" title "{/Symbol f} = 0.55", \
     "data/gamma_phi_0.60.dat" using 1:2 with lines lw 2 lc rgb "#7570b3" title "{/Symbol f} = 0.60", \
     "data/gamma_phi_0.64.dat" using 1:2 with lines lw 2.2 lc rgb "#e7298a" dt 1 title "{/Symbol f} = 0.64 (Warm-Start Only)"

set terminal pngcairo enhanced color font "Helvetica,11" size 1000,700
set output "plots/fig_gamma_profiles.png"
replot
