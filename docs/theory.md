# Theoretical Foundations - OZE_c_solver

This document outlines the statistical mechanics, integral equations, rotational invariant expansions, and numerical methods implemented in `OZE_c_solver`.

---

## 1. Isotropic Ornstein-Zernike (OZ) Equation

For a homogeneous and isotropic fluid, the Ornstein-Zernike integral equation connects the total correlation function $h(r) = g(r) - 1$ to the direct correlation function $c(r)$:

$$h(r) = c(r) + \rho \int c(|\mathbf{r} - \mathbf{r}'|) h(r') d\mathbf{r}'$$

where $\rho = N/V$ is the number density.

In Fourier space, this convolution reduces to an algebraic relation:

$$\hat{h}(k) = \hat{c}(k) + \rho \hat{c}(k) \hat{h}(k) \implies \hat{h}(k) = \frac{\hat{c}(k)}{1 - \rho \hat{c}(k)}$$

The static structure factor $S(k)$ is defined as:

$$S(k) = 1 + \rho \hat{h}(k) = \frac{1}{1 - \rho \hat{c}(k)}$$

### Standard Isotropic Closures

1. **Hypernetted-Chain (HNC)**:
   $$c(r) = \exp[-\beta U(r) + \gamma(r)] - \gamma(r) - 1$$
   where $\gamma(r) = h(r) - c(r)$ is the indirect correlation function.

2. **Rogers-Young (RY)**:
   $$g(r) = \exp[-\beta U(r)] \left( 1 + \frac{\exp[\gamma(r) f(r)] - 1}{f(r)} \right)$$
   with mixing function $f(r) = 1 - e^{-\alpha r}$.

3. **Percus-Yevick (PY)**:
   $$c(r) = g(r) [1 - \exp(\beta U(r))]$$

---

## 2. Non-Spherical & Dipolar Integral Equations

For anisotropic molecules with orientations $\Omega_1, \Omega_2$, pair functions are expanded in rotational invariants $\Phi^{l_1 l_2 l}(12)$:

$$X(12) = \sum_{l_1 l_2 l} X^{l_1 l_2 l}(r) \Phi^{l_1 l_2 l}(\Omega_1, \Omega_2, \hat{r})$$

For dipolar hard spheres (DHS), the primary projections are $(000)$, $(110)$, and $(112)$, representing:
- $000$: Center-of-mass radial correlations.
- $110$: Dipolar angular alignment $\mathbf{s}_1 \cdot \mathbf{s}_2$.
- $112$: Spatial dipole-dipole anisotropy $3(\mathbf{s}_1 \cdot \hat{r})(\mathbf{s}_2 \cdot \hat{r}) - \mathbf{s}_1 \cdot \mathbf{s}_2$.

### Decoupled $\chi$-Mode Reciprocal Space Representation (Blum, 1972)

Under Hankel-Bessel transformation $\tilde{F}^{l_1 l_2 l}(k) = 4\pi i^l \int_0^\infty r^2 j_l(kr) F^{l_1 l_2 l}(r) dr$, the coupled angular OZ equation diagonalizes into decoupled scalar algebraic systems:

- **Isotropic channel ($\chi = \text{iso}$)**:
  $$\tilde{H}^{000}(k) = \frac{\tilde{C}^{000}(k)}{1 - \rho \tilde{C}^{000}(k)}$$

- **Dipolar Mode $\chi = 0$**:
  $$\tilde{C}^0(k) = \tilde{C}^{110}(k) + 2 \tilde{C}^{112}(k) \implies \tilde{H}^0(k) = \frac{\tilde{C}^0(k)}{1 - \frac{\rho}{3} \tilde{C}^0(k)}$$

- **Dipolar Mode $\chi = 1$**:
  $$\tilde{C}^1(k) = \tilde{C}^{110}(k) - \tilde{C}^{112}(k) \implies \tilde{H}^1(k) = \frac{\tilde{C}^1(k)}{1 + \frac{\rho}{3} \tilde{C}^1(k)}$$

---

## 3. Non-Spherical Closure Relations

1. **MSA (Mean Spherical Approximation)**:
   $$c^{000}(r) = c_{\text{HS}}(r), \quad c^{110}(r) = 0, \quad c^{112}(r) = \frac{\beta \mu^2}{r^3} \quad (r > \sigma)$$

2. **LHNC (Linearized HNC)**:
   Linearizes the angular expansion of $\exp[-\beta u + \eta]$ around the isotropic background.

3. **QHNC (Quadratic HNC)**:
   Retains second-order terms via Wigner 9-$j$ angular product couplings.

4. **RHNC (Reference HNC)**:
   Uses exact radial derivative integration with a Percus-Yevick hard-sphere reference system (see [RHNC.md](RHNC.md)).

---

## 4. Thermodynamic Quantities

- **Virial Pressure ($P_v$)**:
  $$\frac{\beta P}{\rho} = 1 - \frac{2\pi\rho}{3} \int_0^\infty r^3 \frac{dU(r)}{dr} g(r) dr$$

- **Isothermal Compressibility ($\chi_T$)**:
  $$\rho k_B T \chi_T = S(k \to 0)$$

- **Internal Energy ($U_{\text{int}}$)**:
  $$\frac{\beta U_{\text{int}}}{N} = 2\pi\rho \int_0^\infty r^2 U(r) g(r) dr$$

---

## 5. References

- **Hansen, J.-P. & McDonald, I. R.** (2013). *Theory of Simple Liquids*. Academic Press.
- **Blum, L.** (1972). *The Journal of Chemical Physics*, 57(5), 1862.
- **Fries, P. H. & Patey, G. N.** (1985). *The Journal of Chemical Physics*, 82(1), 429.
