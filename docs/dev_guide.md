# Developer Guide - OZE_c_solver

This guide is intended for developers and researchers who wish to extend, maintain, or contribute to the **OZE_c_solver** codebase.

---

## 1. Codebase Architecture

```
OZE_c_solver/
├── src/                          # C source files (.c)
│   ├── main.c                    # Entry point, argument parsing, CLI dispatch
│   ├── solver_dipolar.c          # Dipolar OZ solver, Picard iteration, Anderson mixing
│   ├── closures_nonspherical.c   # MSA, LHNC, QHNC, RHNC non-spherical closures
│   ├── structures_nonspherical.c # Hankel & discrete sine transforms for l=0, 2
│   ├── solver_mode2.c            # Coupled mode-2 rotational invariant solver
│   ├── facdes2Y.c                # Isotropic radial OZ solver interface
│   ├── structures.c              # Isotropic potential definitions & Ng acceleration
│   └── math_aux.c                # Fast Fourier transform utilities
├── include/                      # C header files (.h)
│   ├── facdes2Y.h
│   ├── math_aux.h
│   ├── structures.h
│   └── structures_nonspherical.h
├── reports/                      # Benchmark reproduction scripts, figures, LaTeX reports
├── .github/workflows/ci.yml      # GitHub Actions CI workflow
├── Makefile                      # Build automation & test targets
└── README.md                     # Main repository documentation
```

---

## 2. Execution Flow for Dipolar Solver (`solver_dipolar.c`)

1. **CLI Parsing (`src/main.c`)**:
   Reads CLI flags (`--closure`, `--potential`, `--volfactor`, `--temp`, `--dipole`, `--rmax`, `--ramp`, `--init-sk`).
2. **Solver Initialization (`solver_dipolar.c`)**:
   - Computes grid parameters $\Delta r = r_{\max} / N$ and $\Delta k = \pi / r_{\max}$.
   - Precomputes the isotropic hard-sphere Percus-Yevick reference solution ($c_R(r), h_R(r), \eta_R(r), \frac{\partial \ln g_R}{\partial r}$) if closure is RHNC.
   - Initializes $\gamma^{000}, \gamma^{110}, \gamma^{112}$ (either from cold start, continuation stage, or analytical $S(k)$ table).
3. **Iterative Convergence Loop**:
   - **Closure Evaluation (`closures_nonspherical.c`)**: Computes $c^{000}(r), c^{110}(r), c^{112}(r)$ from current $\gamma(r)$ and potential $u(r)$.
   - **Hankel Transforms (`structures_nonspherical.c`)**: Forward transforms $r \cdot c^{l_1 l_2 l}(r) \to k \cdot \tilde{C}^{l_1 l_2 l}(k)$.
   - **$\chi$-Mode Decoupled OZ Algebra**: Solves the decoupled scalar equations for $\tilde{H}^0(k)$ and $\tilde{H}^1(k)$.
   - **Inverse Hankel Transforms**: Inverse transforms $k \cdot \tilde{\Gamma}^{l_1 l_2 l}(k) \to r \cdot \gamma_{\text{new}}^{l_1 l_2 l}(r)$.
   - **Anderson Acceleration**: Mixes up to 4 past iterates to accelerate convergence and suppress oscillations.
4. **Output Writing**:
   Writes final correlation functions and structure factors to `output/`.

---

## 3. Adding a New Interaction Potential

To add a new potential:
1. Open `src/structures.c`.
2. Locate the `POT()` function.
3. Add a new `case ID:` in the `switch(potentialID)`.
4. Implement `U[i*ncols + k]` and derivative `Up[i*ncols + k]` for virial pressure evaluation.
5. Register the potential name in `display_potential_options()` and `print_potential_help()` in `src/main.c`.

---

## 4. Coding Standards

- **Language Standard**: C99.
- **Language Policy**: All comments, identifiers, documentation, and user-facing messages must be in **US English**.
- **Memory Safety**: Every dynamic allocation with `malloc`, `calloc`, or `gsl_vector_alloc` must have a corresponding `free` or `gsl_vector_free`.
- **Testing**: Ensure any new feature compiles cleanly without warnings (`-Wall`) and passes `make test-all`.
