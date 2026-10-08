# plot_speedup.gp
# Speedup factor scaling across volume fractions

set terminal pdfcairo enhanced color font "Helvetica,11" size 4.5,3.0
set output "plots/fig_speedup.pdf"

set title "Warm-Start Speedup Factor vs. Volume Fraction {/Symbol f}\n(Percus-Yevick Hard Spheres, N_{nodes}=4096)" font ",11"
set xlabel "Volume Fraction {/Symbol f}" offset 0,-0.5
set ylabel "Speedup Factor (t_{cold} / t_{warm})"
set grid lc rgb "#e0e0e0" dt 2

set xrange [0.48:0.66]
set xtics 0.50, 0.05, 0.65
set yrange [1.0:7.0]
set ytics 1.0

# Plot speedup for phi <= 0.60, and annotate critical phi = 0.64
set label 1 "Critical Region ({/Symbol f}=0.64):\nCold-Start Diverges\nWarm-Start Converges ({/Symbol \245} Speedup)" at 0.64, 4.5 right font "Helvetica-Bold,8" textcolor rgb "#b22222"
set arrow 1 from 0.638, 4.2 to 0.64, 5.8 lc rgb "#b22222" lw 1.5

plot "data/timing_summary.dat" every ::0::2 using 1:6 with linespoints pt 7 ps 1.3 lc rgb "#377eb8" lw 2 title "Measured Speedup"

set terminal pngcairo enhanced color font "Helvetica,11" size 900,600
set output "plots/fig_speedup.png"
replot
