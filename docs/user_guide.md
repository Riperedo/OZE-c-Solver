# User Guide - OZE_c_solver

Welcome to the **OZE_c_solver** user guide. This document explains how to build, configure, and execute the solver for both isotropic colloidal models and anisotropic dipolar fluids.

---

## 1. Prerequisites and Installation

### Dependencies
Make sure you have `gcc` and the `GSL` (GNU Scientific Library) development packages installed.

```bash
# Ubuntu / Debian
sudo apt-get install -y gcc libgsl-dev build-essential

# Fedora / RHEL
sudo dnf install -y gcc gsl-devel

# macOS
brew install gcc gsl
```

### Compilation
Build the executable by running `make` in the root repository directory:
```bash
make clean && make
```
This produces the executable `build/facdes_solver`.

---

## 2. Basic Execution Syntax

The general CLI syntax is:
```bash
./build/facdes_solver --closure [CLOSURE] --potential [ID] --volfactor [ETA] --temp [T] --nodes [N] [OPTIONS]
```

### Mandatory Arguments

| Argument | Description | Example |
| :--- | :--- | :--- |
| `--closure` | Closure relation: `HNC`, `RY`, `PY`, `MSA`, `LHNC`, `QHNC`, `RHNC` | `RHNC` |
| `--potential` | Numeric potential ID (1 to 16) | `14` |
| `--volfactor` | Packing fraction $\eta = \frac{\pi}{6}\rho\sigma^3$ | `0.418879` |
| `--temp` | Reduced temperature $T^*$ | `1.0` |
| `--nodes` | Number of radial discretization points $N$ (powers of 2 recommended) | `4096` |
| `--knodes` | Number of Fourier space nodes (for spherical outputs) | `1024` |

### Key Optional Flags

| Flag | Description | Default |
| :--- | :--- | :--- |
| `--dipole` | Reduced dipole moment $\mu^*$ ($\sqrt{\beta\mu^2/\sigma^3}$, required for potential 14 & 15) | `0.0` |
| `--rmax` | Spatial domain boundary $r_{\max}/\sigma$ | `15.0` (dipolar) / `10.0` (spherical) |
| `--ramp` | Enables automatic temperature continuation annealing ramp | Disabled |
| `--temp-start` | Initial high temperature for annealing continuation | $2 T^*$ |
| `--temp-steps` | Number of intermediate continuation steps | `10` |
| `--init-sk` | Analytical or prior structure factor `.dat` file for warm-start | None |

---

## 3. Practical Usage Examples

### Example 1: Dipolar Hard Spheres with RHNC Closure
Solve the Reference Hypernetted-Chain equation for dipolar hard spheres at $\rho^* = 0.80$ ($\eta \approx 0.418879$), $\mu^{*2} = 2.0$:
```bash
./build/facdes_solver --closure RHNC --potential 14 \
                      --volfactor 0.418879 --temp 1.0 \
                      --dipole 1.41421356 --nodes 4096 --rmax 30.0
```

### Example 2: High Dipole State via Continuation Ramping
At strong dipole moments ($\mu^{*2} = 2.75$), direct Picard iteration may fail from cold start. Use the continuation ramp:
```bash
./build/facdes_solver --closure RHNC --potential 14 \
                      --volfactor 0.418879 --temp 1.0 \
                      --dipole 1.6583124 --nodes 4096 --rmax 30.0 \
                      --ramp --temp-start 5.0 --temp-steps 15
```

### Example 3: Spherical Hertzian Soft Spheres
```bash
./build/facdes_solver --closure HNC --potential 13 \
                      --volfactor 0.3 --temp 1.0 \
                      --nodes 4096 --knodes 1024
```

---

## 4. Output Files

All calculation results are written to the `output/` directory:

- **Dipolar Systems (Potential 14)**:
  - `output/sk_dipolar_000.dat`, `sk_dipolar_110.dat`, `sk_dipolar_112.dat`: Structure factors $S^{000}(k), S^{110}(k), S^{112}(k)$.
  - `output/gr_dipolar_000.dat`, `gr_dipolar_110.dat`, `gr_dipolar_112.dat`: Pair distributions $g^{000}(r), h^{110}(r), h^{112}(r)$.
  - `output/cr_dipolar_000.dat`, `cr_dipolar_110.dat`, `cr_dipolar_112.dat`: Direct correlation functions $c^{000}(r), c^{110}(r), c^{112}(r)$.

- **Isotropic Systems**:
  - `output/HNC_SdeK.dat` / `output/RY_SdeK.dat`: Static structure factor $S(k)$.
  - `output/HNC_GdeR.dat` / `output/RY_GdeR.dat`: Radial distribution function $g(r)$.

---

## 5. Automated Verification & Reports

You can run automated test targets via `make`:

```bash
make test-all       # Run spherical and non-spherical test suites
make mc-all         # Run full Monte Carlo benchmark pipeline and compile report
make benchmark-all  # Run Wertheim analytical MSA benchmark pipeline
```
