# plot_runtime.gp
# Wall-clock execution time comparison between cold-start and warm-start

set terminal pdfcairo enhanced color font "Helvetica,11" size 4.5,3.0
set output "plots/fig_runtime.pdf"

set title "Execution Time Comparison: Cold-Start vs Warm-Start\n(Percus-Yevick Hard Spheres, N_{nodes}=4096)" font ",11"
set style data histograms
set style histogram cluster gap 1
set style fill solid 0.8 border -1
set boxwidth 0.85

set xlabel "Volume Fraction {/Symbol f}" offset 0,-0.5
set ylabel "Wall-Clock Execution Time (s)"
set grid y lc rgb "#dddddd" dt 2

set yrange [0:35]
set ytics 5
set xtics ("0.50" 0, "0.55" 1, "0.60" 2, "0.64" 3)

set key top left reverse Left

# Annotations for Phi = 0.64 divergence
set label 1 "DIVERGED\n(Limit-Cycle)" at 3.0, 31 center font "Helvetica-Bold,8" textcolor rgb "#b22222"

plot "data/timing_summary.dat" using 2:xtic(1) title "Cold-Start (100-step ramp)" lc rgb "#d95f02", \
     "" using 4 title "Warm-Start ({/Symbol g}(r) seed)" lc rgb "#1b9e77"

set terminal pngcairo enhanced color font "Helvetica,11" size 900,600
set output "plots/fig_runtime.png"
replot
