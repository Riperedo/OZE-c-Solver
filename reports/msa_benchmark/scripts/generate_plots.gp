# Gnuplot script for high-quality MSA benchmark figures
# Terminal settings: PDF with Cairo renderer and enhanced math fonts

set encoding utf8

# Common styling definitions
set style line 1 lc rgb '#0072B2' lt 1 lw 2.5 pt 7   ps 0.7  # Blue solid (Numerical)
set style line 2 lc rgb '#D55E00' lt 2 lw 2.5 pt 6   ps 0.7  # Orange dashed (Analytical)
set style line 3 lc rgb '#009E73' lt 1 lw 2.5 pt 9   ps 0.7  # Green solid
set style line 4 lc rgb '#CC79A7' lt 2 lw 2.5 pt 8   ps 0.7  # Purple dashed
set style line 5 lc rgb '#F0E442' lt 1 lw 2.5 pt 5   ps 0.7  # Yellow solid
set style line 6 lc rgb '#000000' lt 2 lw 2.5 pt 4   ps 0.7  # Black dashed

set grid lc rgb '#E0E0E0' lt 1 lw 0.8
set tics font "Helvetica,11"
set xlabel font "Helvetica,13"
set ylabel font "Helvetica,13"
set key font "Helvetica,10" box opaque

# ----------------------------------------------------------------------
# FIGURE 1: S000(k) comparison for phi in [0.1, 0.3, 0.5] at T=10.0 and T=1.0
# ----------------------------------------------------------------------
set term pdfcairo enhanced color font "Helvetica,12" size 8in, 4.2in
set output "reports/msa_benchmark/plots/fig1_s000_comparison.pdf"

set multiplot layout 1,2 title "Hard-Core Base Projection S^{000}(k): Numerical vs Analytical MSA" font "Helvetica-Bold,14"

# Left Panel: T = 10.0
set title "Temperature T* = 10.0" font "Helvetica,12"
set xlabel "Wavenumber k {/Symbol s}"
set ylabel "S^{000}(k)"
set xrange [0:18]
set yrange [0:3.2]
set key top right

plot "reports/msa_benchmark/data/ana_phi_0.1_T_10.00.dat" u 1:2 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_10.00.dat" u 1:2 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_10.00.dat" u 1:2 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_10.00.dat" u 1:2 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_10.00.dat" u 1:2 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_10.00.dat" u 1:2 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

# Right Panel: T = 1.0
set title "Temperature T* = 1.0" font "Helvetica,12"
plot "reports/msa_benchmark/data/ana_phi_0.1_T_1.00.dat" u 1:2 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_1.00.dat" u 1:2 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_1.00.dat" u 1:2 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_1.00.dat" u 1:2 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_1.00.dat" u 1:2 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_1.00.dat" u 1:2 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

unset multiplot

# ----------------------------------------------------------------------
# FIGURE 2: Patey Projections S110(k) and S112(k) for T=10.0 and T=1.0
# ----------------------------------------------------------------------
set output "reports/msa_benchmark/plots/fig2_patey_projections.pdf"
set multiplot layout 2,2 title "Patey Angular Projections S^{110}(k) and S^{112}(k)" font "Helvetica-Bold,14"

# Panel (1,1): S110 at T=10.0
set title "S^{110}(k) at T* = 10.0" font "Helvetica,12"
set ylabel "S^{110}(k)"
set xrange [0:18]
set yrange [-0.002:0.004]
set key top right
plot "reports/msa_benchmark/data/ana_phi_0.1_T_10.00.dat" u 1:3 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_10.00.dat" u 1:3 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_10.00.dat" u 1:3 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_10.00.dat" u 1:3 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_10.00.dat" u 1:3 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_10.00.dat" u 1:3 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

# Panel (1,2): S112 at T=10.0
set title "S^{112}(k) at T* = 10.0" font "Helvetica,12"
set ylabel "S^{112}(k)"
set yrange [-0.08:0.02]
plot "reports/msa_benchmark/data/ana_phi_0.1_T_10.00.dat" u 1:4 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_10.00.dat" u 1:4 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_10.00.dat" u 1:4 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_10.00.dat" u 1:4 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_10.00.dat" u 1:4 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_10.00.dat" u 1:4 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

# Panel (2,1): S110 at T=1.0
set title "S^{110}(k) at T* = 1.0" font "Helvetica,12"
set ylabel "S^{110}(k)"
set yrange [-0.04:0.08]
plot "reports/msa_benchmark/data/ana_phi_0.1_T_1.00.dat" u 1:3 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_1.00.dat" u 1:3 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_1.00.dat" u 1:3 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_1.00.dat" u 1:3 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_1.00.dat" u 1:3 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_1.00.dat" u 1:3 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

# Panel (2,2): S112 at T=1.0
set title "S^{112}(k) at T* = 1.0" font "Helvetica,12"
set ylabel "S^{112}(k)"
set yrange [-0.4:0.1]
plot "reports/msa_benchmark/data/ana_phi_0.1_T_1.00.dat" u 1:4 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_1.00.dat" u 1:4 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_1.00.dat" u 1:4 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_1.00.dat" u 1:4 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_1.00.dat" u 1:4 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_1.00.dat" u 1:4 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

unset multiplot

# ----------------------------------------------------------------------
# FIGURE 3: Decoupled Chi-modes S^0(k) and S^1(k)
# ----------------------------------------------------------------------
set output "reports/msa_benchmark/plots/fig3_chi_modes.pdf"
set multiplot layout 1,2 title "Invariant Decoupled Chi-Modes S^0(k) and S^1(k) (T* = 1.0)" font "Helvetica-Bold,14"

# Left: S0 (chi=0 mode)
set title "Mode {/Symbol c} = 0 : S^0(k)" font "Helvetica,12"
set xlabel "Wavenumber k {/Symbol s}"
set ylabel "S^0(k)"
set xrange [0:18]
set yrange [0.4:1.4]
set key bottom right
plot "reports/msa_benchmark/data/ana_phi_0.1_T_1.00.dat" u 1:5 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_1.00.dat" u 1:5 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_1.00.dat" u 1:5 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_1.00.dat" u 1:5 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_1.00.dat" u 1:5 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_1.00.dat" u 1:5 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

# Right: S1 (chi=1 mode)
set title "Mode {/Symbol c} = 1 : S^1(k)" font "Helvetica,12"
set ylabel "S^1(k)"
set yrange [0.8:1.6]
set key top right
plot "reports/msa_benchmark/data/ana_phi_0.1_T_1.00.dat" u 1:6 w l ls 2 title "{/Symbol f}=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.1_T_1.00.dat" u 1:6 every 4 w p ls 1 title "{/Symbol f}=0.1 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.3_T_1.00.dat" u 1:6 w l ls 4 title "{/Symbol f}=0.3 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.3_T_1.00.dat" u 1:6 every 4 w p ls 3 title "{/Symbol f}=0.3 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.5_T_1.00.dat" u 1:6 w l ls 6 title "{/Symbol f}=0.5 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.5_T_1.00.dat" u 1:6 every 4 w p ls 5 title "{/Symbol f}=0.5 (Numerical)"

unset multiplot

# ----------------------------------------------------------------------
# FIGURE 4: Thermal Evolution across Temperatures for fixed phi=0.2
# ----------------------------------------------------------------------
set output "reports/msa_benchmark/plots/fig4_thermal_evolution.pdf"
set multiplot layout 1,2 title "Thermal Evolution of Dipolar MSA Structure Factors ({/Symbol f} = 0.2)" font "Helvetica-Bold,14"

# Left: S112 across temperatures
set title "Dipolar Anisotropy S^{112}(k)" font "Helvetica,12"
set xlabel "Wavenumber k {/Symbol s}"
set ylabel "S^{112}(k)"
set xrange [0:18]
set yrange [-1.2:0.2]
set key bottom right
plot "reports/msa_benchmark/data/ana_phi_0.2_T_10.00.dat" u 1:4 w l ls 2 title "T*=10.0 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.2_T_10.00.dat" u 1:4 every 4 w p ls 1 title "T*=10.0 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.2_T_1.00.dat" u 1:4 w l ls 4 title "T*=1.0 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.2_T_1.00.dat" u 1:4 every 4 w p ls 3 title "T*=1.0 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.2_T_0.10.dat" u 1:4 w l ls 6 title "T*=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.2_T_0.10.dat" u 1:4 every 4 w p ls 5 title "T*=0.1 (Numerical)"

# Right: S1 mode across temperatures
set title "Decoupled Mode S^1(k)" font "Helvetica,12"
set ylabel "S^1(k)"
set yrange [0.7:3.5]
set key top right
plot "reports/msa_benchmark/data/ana_phi_0.2_T_10.00.dat" u 1:6 w l ls 2 title "T*=10.0 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.2_T_10.00.dat" u 1:6 every 4 w p ls 1 title "T*=10.0 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.2_T_1.00.dat" u 1:6 w l ls 4 title "T*=1.0 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.2_T_1.00.dat" u 1:6 every 4 w p ls 3 title "T*=1.0 (Numerical)", \
     "reports/msa_benchmark/data/ana_phi_0.2_T_0.10.dat" u 1:6 w l ls 6 title "T*=0.1 (Analytic)", \
     "reports/msa_benchmark/data/num_phi_0.2_T_0.10.dat" u 1:6 every 4 w p ls 5 title "T*=0.1 (Numerical)"

unset multiplot
