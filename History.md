# Development History

## Repository Synchronization and Setup (2026-09-02)

- **Repository Status Verification**:
  - Checked local branch `main` against remote `origin/main` at `https://github.com/Riperedo/OZE-c-Solver.git`.
  - Confirmed local repository is completely up to date with the latest commit (`852f009` - *Enable PY closure in CLI and add wrapper functions*).
  - Validated that remote changes fetch cleanly without discrepancies.

- **Authentication & Credentials Verification**:
  - Validated Git user configuration (`Riperedo` / `rperedo@if.uaslp.mx`).
  - Validated GitHub SSH authentication (`ssh -T git@github.com` successfully authenticated as `Riperedo`).
  - Tested push access via `git push --dry-run origin main`, confirming stored credentials/permissions are active and operational.

- **Untracked Working Tree Items**:
  - Identified untracked script `examples/run_repulsive_Yukawa.sh`.

## Technical Improvement Report Evaluation & Adjustments (2026-09-02)

- Evaluated `docs/reporte_mejoras_oze_solver.md`.
- **Sign Convention Analysis**:
  - Re-examined the sign in `denom1 = 1.0 - (rho / 3.0) * C1;` in `src/solver_dipolar.c`.
  - Analytical validation confirmed that under the basis convention $C^1 = C^{110} - C^{112}$ (where $C^{112} < 0$), the original minus sign is mathematically exact ($1 - (\rho/3)C^1$).
  - Reverted the proposed sign flip, restoring exact parity with the Wertheim/Blum analytical datasets.

- **Implemented Validated Improvements**:
  1. **Dynamic `--rmax` Parameter Support**: Added CLI argument `--rmax <double>` with default $r_{\max} = 15.0$ in `src/main.c`.
  2. **Orthogonal Hankel Transform Quadrature**: Updated `HT2_Direct` and `IHT2_Direct` in `src/math_aux.c` to use uniform quadrature weights `dr` and `dk`, establishing machine-precision discrete adjointness ($< 10^{-13}$).
  3. **Core Gibbs Artifact Elimination**: Precomputed $j_0$ and $j_2$ forward and inverse kernel transformation matrices in `src/solver_dipolar.c`, eliminating high-frequency oscillations in $h^{112}(r)$ near the hard core.
  4. **Anderson Mixing Acceleration**: Implemented multi-vector Anderson acceleration (depth $M=4$) with a small-system linear solver, reducing iteration counts from $\sim 2000$ to $\sim 15 - 38$ iterations.

## Analytical MSA Benchmark & Academic LaTeX Report (2026-09-02)

- **Comprehensive 20-State Thermodynamic Validation**:
  - Benchmarked the numerical solver against all 20 exact analytical MSA structure factor datasets in `test/data_MSA_analitica_HD_dipolar/` at $\mu = 1.0$.
  - State points evaluated:
    - Temperatures: $T^* \in \{10.0, 1.0, 0.1, 0.01\}$
    - Packing fractions: $\phi \in \{0.1, 0.2, 0.3, 0.4, 0.5\}$
- **Micro-Precision Numerical Agreement**:
  - High-temperature regime ($T^* = 10.0$):
    - $\text{RMSE}_{S^{110}} = 2.81 \times 10^{-6}$ ($\phi=0.1$) to $6.32 \times 10^{-5}$ ($\phi=0.5$).
    - $\text{RMSE}_{S^{112}} = 1.05 \times 10^{-4}$ ($\phi=0.1$) to $4.97 \times 10^{-4}$ ($\phi=0.5$).
    - $\text{RMSE}_{S^0} = 2.10 \times 10^{-4}$ ($\phi=0.1$) to $9.94 \times 10^{-4}$ ($\phi=0.5$).
    - $\text{RMSE}_{S^1} = 1.05 \times 10^{-4}$ ($\phi=0.1$) to $5.04 \times 10^{-4}$ ($\phi=0.5$).
  - Moderate coupling regime ($T^* = 1.0$):
    - $\text{RMSE}_{S^{110}} = 2.21 \times 10^{-4}$ to $2.66 \times 10^{-3}$.
    - $\text{RMSE}_{S^{112}} = 9.34 \times 10^{-4}$ to $3.32 \times 10^{-3}$.
    - $\text{RMSE}_{S^0} = 1.88 \times 10^{-3}$ to $7.30 \times 10^{-3}$.
    - $\text{RMSE}_{S^1} = 9.58 \times 10^{-4}$ to $4.13 \times 10^{-3}$.
- **Publication Figures Generated via Gnuplot (Vector PDF)**:
  - `reports/msa_benchmark/plots/fig1_s000_comparison.pdf`: Base hard-core structure factor $S^{000}(k)$.
  - `reports/msa_benchmark/plots/fig2_patey_projections.pdf`: Patey projections $S^{110}(k)$ and $S^{112}(k)$.
  - `reports/msa_benchmark/plots/fig3_chi_modes.pdf`: Decoupled invariant modes $S^0(k)$ and $S^1(k)$.
  - `reports/msa_benchmark/plots/fig4_thermal_evolution.pdf`: Thermal progression from weak to strong dipolar coupling.
  - `reports/msa_benchmark/plots/fig5_error_scaling.pdf`: Error convergence vs. packing fraction $\phi$.
- **Formal Academic Report in LaTeX**:
  - Generated and compiled [`reports/msa_benchmark/msa_benchmark_report.pdf`](file:///home/jinzo/Desktop/codigos/OZE_c_solver/reports/msa_benchmark/msa_benchmark_report.pdf) (9 pages, fully structured with mathematical formulation, numerical algorithms, error metric tables, embedded figures, and physical discussion).

## Temperature Continuation (Annealing) & Robust Low-Temperature Solver (2026-09-02)

- **Geometric Temperature Continuation Framework**:
  - Implemented multi-stage temperature continuation schedule:
    $$T_s^* = T_{\text{start}}^* \left( \frac{T_{\text{target}}^*}{T_{\text{start}}^*} \right)^{\frac{s}{N-1}}, \quad s = 0, \dots, N-1$$
  - Seamless warm-starting transferring converged correlation matrices $\mathbf{c}_s(r)$ and indirect screening profiles $\bm{\eta}_s(r)$ across temperature stages.
  - Multi-tier adaptive tolerance schedule ($\varepsilon = 10^{-4}$ for intermediate stages, $\varepsilon = 10^{-6}$ for target temperature).
- **Scale-Invariant Regularized Anderson Acceleration**:
  - Upgraded Anderson mixing to the unconstrained difference-vector Walker--Ni formulation.
  - Added scale-invariant Tikhonov regularization $\lambda = 10^{-4} \max_j M_{jj}$ to the $(m-1) \times (m-1)$ normal equation matrix.
  - Added smooth weight bounding to prevent over-extrapolation while avoiding fallback stalls to unaccelerated Picard steps.
  - Reduced per-stage iteration counts to $\sim 8 - 14$ iterations, traversing from $T^* = 10.0$ down to $T^* = 0.10$ ($\beta\mu^2 = 10.0$) in under 3 seconds total.
- **CLI Options Added to `src/main.c`**:
  - `--temp-start <double>`: Initial temperature for continuation path (e.g., `10.0`).
  - `--temp-steps <int>`: Number of logarithmic stages (default `10`).
  - `--ramp`: Activates automatic continuation schedule with smart high-temperature default.
- **Academic Report Updated & Synchronized**:
  - Added Section 4.5 and Table II to [`reports/msa_benchmark/msa_benchmark_report.tex`](file:///home/jinzo/Desktop/codigos/OZE_c_solver/reports/msa_benchmark/msa_benchmark_report.tex) detailing the temperature continuation algorithm and convergence characteristics.
  - Re-ran full benchmark suite across all 20 states using temperature continuation, regenerating benchmark datasets in `reports/msa_benchmark/data/` and vector figures.
  - Characterized finite-box $R_{\max}$ dipolar truncation effects ($1/r^3$ tail discretization) at strong coupling ($\beta\mu^2 \ge 10.0$).
  - Recompiled LaTeX report to [`reports/msa_benchmark/msa_benchmark_report.pdf`](file:///home/jinzo/Desktop/codigos/OZE_c_solver/reports/msa_benchmark/msa_benchmark_report.pdf).


