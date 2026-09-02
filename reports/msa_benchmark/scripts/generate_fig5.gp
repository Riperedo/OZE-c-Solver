# Gnuplot script for error scaling using named datablocks
set encoding utf8
set term pdfcairo enhanced color font "Helvetica,12" size 8in, 4.2in
set output "reports/msa_benchmark/plots/fig5_error_scaling.pdf"

$ERR_T10 << EOD
# phi RMSE_S000 RMSE_S110 RMSE_S112 RMSE_S10 RMSE_S11
0.1 4.27e-03 2.81e-06 1.05e-04 2.10e-04 1.05e-04
0.2 7.79e-03 1.10e-05 2.07e-04 4.13e-04 2.08e-04
0.3 1.39e-02 2.41e-05 3.06e-04 6.12e-04 3.09e-04
0.4 2.67e-02 4.12e-05 4.03e-04 8.06e-04 4.07e-04
0.5 5.54e-02 6.32e-05 4.97e-04 9.94e-04 5.04e-04
EOD

$ERR_T1 << EOD
# phi RMSE_S000 RMSE_S110 RMSE_S112 RMSE_S10 RMSE_S11
0.1 4.27e-03 2.21e-04 9.34e-04 1.88e-03 9.58e-04
0.2 7.79e-03 6.97e-04 1.67e-03 3.45e-03 1.77e-03
0.3 1.39e-02 1.29e-03 2.28e-03 4.84e-03 2.52e-03
0.4 2.67e-02 1.94e-03 2.82e-03 6.11e-03 3.29e-03
0.5 5.54e-02 2.66e-03 3.32e-03 7.30e-03 4.13e-03
EOD

set multiplot layout 1,2 title "Numerical Error Convergence vs Packing Fraction {/Symbol f}" font "Helvetica-Bold,14"

set grid lc rgb '#E0E0E0' lt 1 lw 0.8
set tics font "Helvetica,11"
set xlabel "Packing Fraction {/Symbol f}" font "Helvetica,13"
set ylabel "Root Mean Square Error (RMSE)" font "Helvetica,13"
set logscale y
set format y "10^{%T}"
set key top left font "Helvetica,10" box opaque
set xrange [0.08:0.52]

# Left: S000 and S10 errors vs phi for T=10.0 and T=1.0
set title "Hard-Core & Isotropic Mode Errors" font "Helvetica,12"
plot $ERR_T10 u 1:2 w lp lc rgb '#0072B2' pt 7 lw 2 ps 1 title "S^{000} (T*=10.0)", \
     $ERR_T1  u 1:2 w lp lc rgb '#D55E00' pt 9 lw 2 ps 1 title "S^{000} (T*=1.0)", \
     $ERR_T10 u 1:5 w lp lc rgb '#009E73' pt 5 lw 2 ps 1 title "S^{0} (T*=10.0)", \
     $ERR_T1  u 1:5 w lp lc rgb '#CC79A7' pt 11 lw 2 ps 1 title "S^{0} (T*=1.0)"

# Right: S110 and S112 errors vs phi for T=10.0 and T=1.0
set title "Anisotropic Patey Mode Errors" font "Helvetica,12"
plot $ERR_T10 u 1:3 w lp lc rgb '#0072B2' pt 7 lw 2 ps 1 title "S^{110} (T*=10.0)", \
     $ERR_T1  u 1:3 w lp lc rgb '#D55E00' pt 9 lw 2 ps 1 title "S^{110} (T*=1.0)", \
     $ERR_T10 u 1:4 w lp lc rgb '#009E73' pt 5 lw 2 ps 1 title "S^{112} (T*=10.0)", \
     $ERR_T1  u 1:4 w lp lc rgb '#CC79A7' pt 11 lw 2 ps 1 title "S^{112} (T*=1.0)"

unset multiplot
