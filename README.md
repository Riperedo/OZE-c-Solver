# OZE_c_solver: Ornstein-Zernike Equation Solver

[![CI](https://github.com/Riperedo/OZE-c-Solver/actions/workflows/ci.yml/badge.svg)](https://github.com/Riperedo/OZE-c-Solver/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Language: C](https://img.shields.io/badge/Language-C99-orange.svg)](https://en.wikipedia.org/wiki/C99)

A high-performance C solver for the **Ornstein-Zernike (OZ) integral equation** in classical and colloidal statistical mechanics. It computes the total correlation functions $h(r)$, direct correlation functions $c(r)$, radial distribution functions $g(r)$, and static structure factors $S(k)$.

The package supports standard isotropic simple fluids as well as generalized **anisotropic (non-spherical / dipolar)** fluids via rotational invariant projections $\chi$-space decoupled linear systems ($000, 110, 112$).

---

## 🌟 Key Features

- **Isotropic Closures**:
  - **HNC** (Hypernetted-Chain)
  - **RY** (Rogers-Young with thermodynamic self-consistency)
  - **PY** (Percus-Yevick)
- **Non-Spherical (Dipolar Hard Spheres) Closures**:
  - **MSA** (Mean Spherical Approximation, Blum & Wertheim formulation)
  - **LHNC** (Linearized Hypernetted-Chain)
  - **QHNC** (Quadratic Hypernetted-Chain)
  - **RHNC** (Reference Hypernetted-Chain with hard-sphere reference system and exact Fries & Patey radial derivatives)
- **Numerical Acceleration & Robustness**:
  - Ng acceleration and multi-dimensional Anderson mixing ($M=4$)
  - Geometric temperature continuation annealing (`--ramp`)
  - Analytical $S(k)$ warm-starting (`--init-sk`)
  - Hankel transform kernels with fast sine transforms ($O(N \log N)$)
- **Academic Benchmark Pipelines**:
  - Full automated validation against Wertheim's exact analytical dipolar solution
  - Full automated validation against Monte Carlo simulation benchmarks (Fries & Patey 1985, Patey-Levesque-Weis 1979)
  - Publication-ready LaTeX report generators (`make report-mc`, `make report`)

---

## 📋 Requirements & Dependencies

- **C Compiler**: GCC (recommended $\ge 9.0$) or Clang with C99 support
- **Scientific Libraries**: [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
- **Python / Plotting** (Optional, for benchmark suites and reports): Python 3.8+, NumPy, Matplotlib, Gnuplot, pdfLaTeX

### Installing Dependencies

```bash
# Ubuntu / Debian
sudo apt-get update
sudo apt-get install -y gcc libgsl-dev python3 python3-numpy python3-matplotlib gnuplot texlive-latex-base texlive-latex-extra

# Fedora / RHEL
sudo dnf install -y gcc gsl-devel python3 python3-numpy python3-matplotlib gnuplot texlive-scheme-basic

# macOS (Homebrew)
brew install gcc gsl python gnuplot
```

---

## 🚀 Quick Start

### 1. Compilation
```bash
make clean
make
```

### 2. Running Verification Tests
```bash
make test-all       # Runs both spherical and non-spherical test suites
```

---

## 💻 Command-Line Interface (CLI)

```bash
./build/facdes_solver --closure <CLOSURE> --potential <ID> [OPTIONS]
```

### Core CLI Arguments

| Argument | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `--closure` | `string` | Thermodynamic closure (`HNC`, `RY`, `PY`, `MSA`, `LHNC`, `QHNC`, `RHNC`) | `RHNC` |
| `--potential` | `int` | Interaction potential ID (1 to 16) | `14` |
| `--volfactor` | `double` | Packing fraction $\eta = \frac{\pi}{6} \rho \sigma^3$ | `0.418879` |
| `--temp` | `double` | Reduced temperature $T^*$ (or $1/\beta$) | `1.0` |
| `--nodes` | `int` | Real-space spatial grid resolution $N$ | `4096` |
| `--knodes` | `int` | Number of Fourier space output nodes (optional for dipolar) | `1024` |

### Optional & Advanced Flags

| Flag | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `--dipole` | `double` | Reduced dipole moment $\mu^*$ ($\mu^{*2} = \beta \mu^2 / \sigma^3$) | `0.0` |
| `--rmax` | `double` | Radial cutoff domain $r_{\max} / \sigma$ | `15.0` (dipolar) / `10.0` (spherical) |
| `--ramp` | flag | Enable automated geometric temperature continuation ramp | Disabled |
| `--temp-start` | `double` | High initial temperature for annealing continuation | $2 T^*$ |
| `--temp-steps` | `int` | Number of intermediate continuation temperature stages | `10` |
| `--init-sk` | `string` | Input `.dat` file with prior structure factors for warm-start | None |
| `--lambda_a` | `double` | Attractive / potential range parameter | `0.0` |
| `--lambda_r` | `double` | Repulsive / screening parameter | `0.0` |
| `--temp2` | `double` | Secondary temperature / step width | `1.0` |

---

## 🧪 Interaction Potentials Catalog

| ID | Potential Name | Functional Form / Physics |
| :-: | :--- | :--- |
| **1** | Inverse Power Law (IPL) | $U(r) = T^* (\sigma / r)^\lambda$ |
| **2** | Truncated Lennard-Jones | Repulsive LJ truncated at potential minimum |
| **3** | Truncated Lennard-Jones 2 | LJ potential with minimum shifted to $r=\sigma$ |
| **4** | Double Yukawa | Attractive + Repulsive Yukawa fluid |
| **5** | Attractive Yukawa | $U(r) = -T^* \frac{\exp[-\lambda (r-1)]}{r}$ |
| **6** | Repulsive Yukawa | $U(r) = T^* \frac{\exp[-\lambda (r-1)]}{r}$ |
| **7** | Hard Sphere (HS) | $U(r) = \infty$ ($r < \sigma$), $0$ ($r \ge \sigma$) |
| **8** | Square Shoulder | Step function barrier for $\sigma < r < T_2$ |
| **9** | Down-Hill Potential | Linearly decreasing repulsive shoulder |
| **10** | Gaussian Core Model | $U(r) = T^* \exp(-(r/\sigma)^2)$ |
| **11** | Linear Ramp | Linearly decreasing step |
| **12** | Soft Core Step | $U(r) = E (1 - r/\sigma)^n$ |
| **13** | Hertzian Soft Spheres | $U(r) = E (1 - r/\sigma)^{5/2} \Theta(\sigma - r)$ |
| **14** | **Dipolar Hard Spheres (DHS)** | Hard cores + point dipoles ($000, 110, 112$ projections) |
| **15** | DHS Extended (Modes 2) | DHS coupled up to order $m, n \le 2$ |
| **16** | Soft Shoulder | Smooth hyperbolic tangent shoulder |

---

## 📖 Examples

### 1. Dipolar Hard Spheres with RHNC Closure (Patey State)
Solve the Reference Hypernetted-Chain equation for dipolar hard spheres at $\rho^* = 0.80$ ($\eta \approx 0.418879$), $\mu^{*2} = 2.0$:
```bash
./build/facdes_solver --closure RHNC --potential 14 \
                      --volfactor 0.418879 --temp 1.0 \
                      --dipole 1.41421356 --nodes 4096 --rmax 30.0
```

### 2. Dipolar DHS with Temperature Continuation Ramping
Solve high-dipole regime ($\mu^{*2} = 2.75$) via annealing from $T^*=5.0$:
```bash
./build/facdes_solver --closure RHNC --potential 14 \
                      --volfactor 0.418879 --temp 1.0 \
                      --dipole 1.6583124 --nodes 4096 --rmax 30.0 \
                      --ramp --temp-start 5.0 --temp-steps 15
```

### 3. Spherical Fluid (Hertzian Model)
```bash
./build/facdes_solver --closure HNC --potential 13 \
                      --volfactor 0.3 --temp 1.0 \
                      --nodes 4096 --knodes 1024
```

---

## 📊 Output Files

Outputs are saved in the `output/` directory:

- **Isotropic fluids**:
  - `output/HNC_SdeK.dat`: Static structure factor $S(k)$.
  - `output/HNC_GdeR.dat`: Radial distribution function $g(r)$.
- **Dipolar fluids (Potential 14)**:
  - `output/sk_dipolar_000.dat`, `sk_dipolar_110.dat`, `sk_dipolar_112.dat`: Projections $S^{000}(k)$, $S^{110}(k)$, and $S^{112}(k)$.
  - `output/gr_dipolar_000.dat`, `gr_dipolar_110.dat`, `gr_dipolar_112.dat`: Projections $g^{000}(r)$, $h^{110}(r)$, and $h^{112}(r)$.
  - `output/cr_dipolar_000.dat`, `cr_dipolar_110.dat`, `cr_dipolar_112.dat`: Direct correlation projections $c^{l_1 l_2 l}(r)$.

---

## 🔬 Benchmark Reproducibility & Reports

The repository includes complete automated reproduction pipelines:

### Percus-Yevick Hard Sphere Analytical Benchmark Suite
Reproduces the exact Wertheim-Thiele analytical comparison across all 7 packing fractions ($\phi \in [0.10, 0.60]$):
```bash
make py-all
```
This runs the full solver suite, generates publication figures, and compiles [`reports/py_hard_sphere_benchmark/py_hard_sphere_report.pdf`](reports/py_hard_sphere_benchmark/py_hard_sphere_report.pdf).

### Monte Carlo Simulation Benchmark Suite
Reproduces Fries & Patey (1985) and Patey, Levesque & Weis (1979) simulations across 4 density regimes ($\rho^* \in [0.15, 0.80]$):
```bash
make mc-all
```
This generates all comparative datasets, publication figures, and compiles [`reports/monte_carlo_benchmark/monte_carlo_benchmark_report.pdf`](reports/monte_carlo_benchmark/monte_carlo_benchmark_report.pdf).

### Wertheim Analytical MSA Benchmark Suite
```bash
make benchmark-all
make report
```
Compiles [`reports/msa_benchmark/msa_benchmark_report.pdf`](reports/msa_benchmark/msa_benchmark_report.pdf).

---

## 📁 Repository Structure

```
OZE_c_solver/
├── src/                          # C99 source code
│   ├── main.c                    # CLI entry point and argument parsing
│   ├── solver_dipolar.c          # Dipolar OZ solver (MSA, LHNC, QHNC, RHNC)
│   ├── closures_nonspherical.c   # Non-spherical closure equations & PY reference
│   ├── structures_nonspherical.c # Hankel & discrete sine transforms
│   ├── solver_mode2.c            # Coupled mode-2 extension solver
│   ├── facdes2Y.c                # Isotropic radial OZ solver core
│   ├── math_aux.c                # Mathematical utilities & Ng acceleration
│   └── structures.c              # Isotropic potentials and thermodynamic routines
├── include/                      # Header declarations
│   ├── facdes2Y.h
│   ├── math_aux.h
│   ├── structures.h
│   └── structures_nonspherical.h
├── reports/                      # Academic benchmark reports & reproduction scripts
│   ├── py_hard_sphere_benchmark/ # Percus-Yevick exact analytical report, plots, data
│   ├── monte_carlo_benchmark/    # MC comparison report, scripts, plots, and tables
│   └── msa_benchmark/            # Wertheim analytical MSA report and scripts
├── test/                         # Reference datasets (Monte Carlo & analytical)
│   └── data_gdr_PY_analitica/    # Exact Percus-Yevick HS datasets (phi=0.10 to 0.60)
├── data/                         # Analytical S(k) tables for warm-start
├── docs/                         # Extended mathematical & technical documentation
├── .github/workflows/ci.yml      # Continuous integration workflow
├── Makefile                      # Build automation & benchmark targets
├── LICENSE                       # MIT License
└── README.md                     # Project documentation
```

---

## 📚 References

1. **Fries, P. H. & Patey, G. N.** (1985). *The solution of the RHNC equation for dipolar hard spheres and the calculation of the dielectric constant*. The Journal of Chemical Physics, 82(1), 429–446.
2. **Patey, G. N., Levesque, D., & Weis, J. J.** (1979). *Structure of dipolar fluids: A comparison of theoretical and Monte Carlo results*. Molecular Physics, 38(1), 219–232.
3. **Wertheim, M. S.** (1971). *Exact solution of the mean spherical model for fluids of hard spheres with dipoles*. Physical Review Letters, 27(26), 1776.
4. **Blum, L.** (1972). *Mean spherical model for asymmetric electrolytes and dipolar hard spheres*. The Journal of Chemical Physics, 57(5), 1862–1869.
5. **Hansen, J.-P. & McDonald, I. R.** (2013). *Theory of Simple Liquids* (4th ed.). Academic Press.

---

## 👥 Authors & Acknowledgments

- **Ricardo Peredo Ortiz** — *Institute of Physics, Autonomous University of San Luis Potosí (UASLP)*
- **Jonathan Josué Elisea Espinoza** — *Institute of Physics, Autonomous University of San Luis Potosí (UASLP)*
- **Non-Equilibrium Soft Matter Group** — *Instituto de Física, UASLP*

---

## 📄 License

This project is open-source software licensed under the [MIT License](LICENSE).
