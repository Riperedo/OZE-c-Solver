---
name: Ornstein-Zernike Non-Spherical Solver
description: Domain knowledge, numerical best practices, and scripts for extending the Ornstein-Zernike equation solver for generalized non-spherical potentials (dipoles, quadrupoles, etc.).
---

# Ornstein-Zernike (OZE) Non-Spherical Solver Skill

This skill provides the domain-specific knowledge and best practices required to develop, extend, and debug the Ornstein-Zernike solver for non-spherical particles using rotational invariant projections.

## 1. Domain Knowledge & Notations

When working with non-spherical OZE closures (MSA, LHNC, QHNC, RHNC), two coordinate frameworks are primarily mixed in the literature, which can cause severe numerical divergences if not carefully converted:

*   **Fries & Patey (F&P) Convention**: The standard basis used to output real-space radial projections like $h^{000}(r), h^{110}(r), h^{112}(r)$. These are the physical components mapped to the coordinate frame of the interacting particles.
*   **Blum's Trace Formalism**: Used to algebraically decouple the OZE in $k$-space into independent matrices for each $\chi$-mode. 

### Critical Normalization Trap ($P$ Matrix)
To decouple the matrices, you must map from F&P to Blum and back. The matrix inverse operation for a $\chi$-mode requires a coupling matrix $P$:
$N^\chi = (-1)^\chi \rho C^\chi P (I - (-1)^\chi \rho C^\chi P)^{-1} C^\chi$

*   The components of the F&P vector $\tilde{C}^{m n, \chi}$ must be divided by an exact normalization factor $y^{mnl} = \sqrt{2m+1}\sqrt{2n+1} \begin{pmatrix} m & n & l \\ 0 & 0 & 0 \end{pmatrix}$ before entering Blum's trace matrix.
*   The $P_{\mu\nu}$ matrix evaluates to $(-1)^\mu \delta_{\mu, -\nu}$. Since the equations are restricted to a single $\chi$ where $\mu = \nu = \chi$, $P$ is strictly diagonal but carries alternating signs depending on parity.

## 2. Numerical Integration (Hankel / Bessel Transforms)

Direct integration of $j_l(kr)$ over $N^2$ loops using standard integration (e.g., Simpson's or Trapezoidal) will accumulate $O(1)$ orthogonality drift for $l > 0$, guaranteeing that the Picard iteration will explode to infinity.

**Best Practice:**
Do **NOT** use numerical integration for $j_l(kr)$. Instead, strictly map the transforms into Discrete Sine Transforms (DST) over linear grids ($r_j = j \Delta r, k_i = i \Delta k$):
*   **l = 0:** Exactly maps to $\sin(kr) / kr$.
*   **l > 0:** Must be analytically expanded using precise derivatives of $\sin$ and $\cos$. Use explicit algebraic expansions like:
    *   $l=1$: $(\sin(kr)/kr^2) - (\cos(kr)/kr)$
    *   $l=2$: $((3/kr^3) - (1/kr))\sin(kr) - (3/kr^2)\cos(kr)$

## 3. Python Code Generators (`aux/` directory)

Because the Wigner 3-j algebra and symbolic inverse matrices are extremely error-prone to write by hand in C, always use Python generators `sympy.physics.wigner` to write the C functions.

The following reference scripts are available in the project's `aux/` directory:
*   `generate_mode2_solver_p.py`: Generates the exact $C$ code for decoupling the 14 even-parity modes (Potential 15) and analytically inverting the $P$-coupled matrices.
*   `generate_mode2_transforms.py`: Generates the exact algebraic C expansions for the required spherical Bessel functions $j_l(kr)$ up to $l=4$.

## 4. Iteration Stability (Picard Method)

*   **Mixing Parameter ($\alpha$)**: Begin with strongly damped Picard iterations ($\alpha \le 0.3$) for non-spherical modes.
*   **Error Metric**: Track the L2 error across all projections simultaneously ($\sqrt{\sum (c_{new} - c_{old})^2 / (N_{proj} \times N_{nodes})}$). If the error jumps to order $10^5$ within 50 iterations, it indicates either a missing $P$ matrix (-1) factor or orthogonality drift in your Hankel transform.
