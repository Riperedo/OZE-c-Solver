# Reference Hypernetted-Chain (RHNC) Formulation for Dipolar Hard Spheres

This document details the mathematical theory, rotational invariant projections, and exact radial derivative formulation of the Reference Hypernetted-Chain (RHNC) closure for dipolar hard-sphere fluids, following the classic formulation by Fries & Patey (1985).

---

## 1. Perturbation Difference Variables

The reference fluid theory of Lado (1973), adapted for anisotropic dipolar hard spheres by Fries & Patey (1985), decomposes all correlation functions into an isotropic reference contribution (subscript $R$, chosen as the pure Percus-Yevick hard-sphere fluid) and an anisotropic perturbation (difference operator $\Delta$):

$$\Delta X(12) = X(12) - X_R(r)$$

Because the non-spherical angular projections ($110, 112$) vanish identically in the isotropic hard-sphere reference system ($c_R^{mnl} = 0, h_R^{mnl} = 0, \eta_R^{mnl} = 0$ for $mnl \neq 000$), the difference variables satisfy:

- **For Anisotropic Channels ($mnl \neq 000$):**
  $$\Delta c^{mnl}(r) = c^{mnl}(r)$$
  $$\Delta h^{mnl}(r) = h^{mnl}(r)$$
  $$\Delta \eta^{mnl}(r) = \eta^{mnl}(r)$$
  $$\Delta u^{mnl}(r) = u^{mnl}(r) = -\frac{\mu^2}{r^3} \quad (\text{for channel } 112)$$

- **For the Isotropic Channel ($000$):**
  $$\Delta c^{000}(r) = c^{000}(r) - c_R(r)$$
  $$\Delta h^{000}(r) = h^{000}(r) - h_R(r)$$
  $$\Delta \eta^{000}(r) = \eta^{000}(r) - \eta_R(r)$$
  $$\Delta u^{000}(r) = u^{000}(r) - u_R(r) = 0$$

---

## 2. Core Boundary Condition ($r < \sigma$)

Inside the hard core, particle impenetrability enforces $g(12) = 0$ and $g_R(r) = 0$, meaning $h(12) = -1$ and $h_R(r) = -1$. Consequently, $\Delta h(12) = 0$ inside the core. The Ornstein-Zernike equation reduces to:

$$\Delta c^{mnl}(r) = -\Delta \eta^{mnl}(r) \quad (r < \sigma)$$

This exact linear relation applies to all channels ($000, 110, 112$) inside the core.

---

## 3. Exact Radial Derivative RHNC Closure ($r > \sigma$)

Direct angular projection of the exponential RHNC bridge term is analytically intractable. Fries & Patey (1985) differentiate the RHNC closure with respect to $r$:

$$\frac{\partial \Delta c(12)}{\partial r} = -\Delta h(12) \frac{\partial \Delta W(12)}{\partial r} - h_R(r) \frac{\partial \Delta W(12)}{\partial r} + \Delta h(12) \frac{\partial \ln g_R(r)}{\partial r} - \beta \frac{\partial \Delta u(12)}{\partial r}$$

where the perturbation potential of mean force is:
$$\Delta W(12) = -\Delta \eta(12) + \beta \Delta u(12)$$

Projecting onto the rotational invariant basis $\Phi^{000}, \Phi^{110}, \Phi^{112}$ yields the exact system of ordinary differential equations:

### Isotropic Mode ($000$):
$$\left( \frac{\partial \Delta c}{\partial r} \right)^{000} = h^{000}(r) X'^{000}(r) + \frac{1}{3} h^{110}(r) X'^{110}(r) + \frac{2}{3} h^{112}(r) X'^{112}(r) + [h^{000}(r) - h_R(r)] \frac{\partial \ln g_R(r)}{\partial r}$$

### Dipolar Invariant Mode ($110$):
$$\left( \frac{\partial \Delta c}{\partial r} \right)^{110} = h^{000}(r) X'^{110}(r) + h^{110}(r) X'^{000}(r) + h^{110}(r) \frac{\partial \ln g_R(r)}{\partial r}$$

### Dipolar Invariant Mode ($112$):
$$\left( \frac{\partial \Delta c}{\partial r} \right)^{112} = h^{000}(r) X'^{112}(r) + h^{112}(r) X'^{000}(r) + h^{112}(r) \frac{\partial \ln g_R(r)}{\partial r} - \frac{3 \beta \mu^2}{r^4}$$

where:
- $X'^{000}(r) = \frac{\partial \eta^{000}}{\partial r} - \frac{\partial \eta_R}{\partial r}$
- $X'^{110}(r) = \frac{\partial \eta^{110}}{\partial r}$
- $X'^{112}(r) = \frac{\partial \eta^{112}}{\partial r} - \frac{3 \beta \mu^2}{r^4}$

---

## 4. Backwards Cumulative Integration

Because $\Delta c^{mnl}(r) \to 0$ as $r \to \infty$, the direct correlation functions are obtained by backwards integration from $r_{\max}$:

$$\Delta c^{mnl}(r) = -\int_r^{r_{\max}} \left( \frac{\partial \Delta c^{mnl}}{\partial s} \right) ds$$

Evaluated via the trapezoidal rule:
$$\Delta c^{mnl}(r_i) = \Delta c^{mnl}(r_{i+1}) + \frac{\Delta r}{2} \left[ \left( \frac{\partial \Delta c^{mnl}}{\partial r} \right)_{r_i} + \left( \frac{\partial \Delta c^{mnl}}{\partial r} \right)_{r_{i+1}} \right]$$

---

## 5. References

- **Fries, P. H. & Patey, G. N.** (1985). *The solution of the RHNC equation for dipolar hard spheres and the calculation of the dielectric constant*. *The Journal of Chemical Physics*, 82(1), 429–446.
- **Lado, F.** (1973). *Integral equations for fluids of hard spheres with long-range interactions*. *Physical Review A*, 8(5), 2548.
