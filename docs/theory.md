# Fundamentos Teóricos - OZE_c_solver

Este documento describe brevemente la física y los métodos numéricos implementados en el solver.

## 1. Ecuación de Ornstein-Zernike (OZ)

La ecuación integral de Ornstein-Zernike es fundamental en la teoría de líquidos para relacionar la función de correlación directa $c(r)$ con la función de correlación total $h(r) = g(r) - 1$.

Para un sistema isotrópico y homogéneo:

$$ h(r) = c(r) + \rho \int c(|\mathbf{r} - \mathbf{r}'|) h(r') d\mathbf{r}' $$

Donde:
- $h(r)$: Función de correlación total.
- $c(r)$: Función de correlación directa.
- $\rho$: Densidad numérica del sistema.

En el espacio de Fourier, esta ecuación se simplifica a una relación algebraica:

$$ \hat{h}(k) = \hat{c}(k) + \rho \hat{c}(k) \hat{h}(k) $$

$$ \hat{h}(k) = \frac{\hat{c}(k)}{1 - \rho \hat{c}(k)} $$

El factor de estructura estático $S(k)$ se relaciona con $\hat{h}(k)$ mediante:

$$ S(k) = 1 + \rho \hat{h}(k) = \frac{1}{1 - \rho \hat{c}(k)} $$

## 2. Relaciones de Cierre

Para resolver la ecuación OZ, que tiene dos incógnitas ($h$ y $c$), se necesita una ecuación adicional llamada "relación de cierre".

### Hypernetted Chain (HNC)
Aproximación útil para sistemas con interacciones de largo alcance (e.g., Coulomb, Yukawa).

$$ c(r) = \exp[-\beta U(r) + \gamma(r)] - \gamma(r) - 1 $$

Donde $\gamma(r) = h(r) - c(r)$ es la función de correlación indirecta.

### Rogers-Young (RY)
Un cierre híbrido que interpola entre Percus-Yevick (PY) y HNC. Es termodinámicamente consistente para esferas duras y potenciales repulsivos suaves.

$$ g(r) = \exp[-\beta U(r)] \left( 1 + \frac{\exp[\gamma(r) f(r)] - 1}{f(r)} \right) $$

Con la función de mezcla:
$$ f(r) = 1 - e^{-\alpha r} $$

El parámetro $\alpha$ se ajusta iterativamente para asegurar la consistencia entre la presión virial y la compresibilidad.

## 3. Método Numérico

El solver utiliza el **método de Ng** para acelerar la convergencia de la solución iterativa.

1.  **Inicialización**: Se comienza con una suposición inicial para $\gamma(r)$ (usualmente 0).
2.  **Iteración de Picard**: Se calcula $c(r)$ usando el cierre, luego se transforma a Fourier para obtener $\hat{c}(k)$, se usa OZ para obtener $\hat{\gamma}(k)$, y se transforma inversamente para obtener un nuevo $\gamma(r)$.
3.  **Aceleración de Ng**: En lugar de usar simplemente el resultado de la última iteración, el método de Ng utiliza una combinación lineal de las iteraciones anteriores para predecir una solución más cercana a la convergencia, minimizando el residuo.

### Transformada de Fourier
Se utiliza la Transformada Rápida de Fourier (FFT) para alternar eficientemente entre el espacio real y el recíproco. Debido a la simetría esférica, el problema se reduce a transformadas seno unidimensionales.

## 4. Propiedades Termodinámicas

Una vez obtenidas las funciones de correlación, se calculan propiedades macroscópicas:

- **Presión Virial ($P_v$)**:
  $$ \frac{\beta P}{\rho} = 1 - \frac{2\pi\rho}{3} \int_0^\infty r^3 \frac{dU(r)}{dr} g(r) dr $$

- **Compresibilidad Isotérmica ($\chi_T$)**:
  $$ \rho k_B T \chi_T = S(k \to 0) $$

- **Energía Interna ($U_{int}$)**:
  $$ \frac{\beta U_{int}}{N} = 2\pi\rho \int_0^\infty r^2 U(r) g(r) dr $$
