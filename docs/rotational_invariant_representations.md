# Representaciones de Invariantes Rotacionales para Fluidos Dipolares

## Convención utilizada en este código

Este código usa la **convención de Patey (1977/1980)**, que expande la función de correlación de par $h(12)$ como:

$$h(12) = h^{000}(r) + h^{110}(r)\,\Phi^{110}(12) + h^{112}(r)\,\Phi^{112}(12)$$

donde los índices $(m,n,l)$ son los órdenes angulares de la expansión en armónicos esféricos generalizados (Wigner).

---

## Comparación de las tres representaciones principales

### 1. Wertheim (1971)

Wertheim expande $h(12)$ en términos de productos escalares geométricos de los vectores unitarios dipolar $\hat{s}_1, \hat{s}_2$ y de la dirección intermolecular $\hat{r}$:

$$h(12) = h_s(r) + h_\Delta(r)\,\Delta(12) + h_D(r)\,D(12)$$

con las funciones angulares:

| Función      | Definición                                                                          |
| ------------ | ----------------------------------------------------------------------------------- |
| $\Delta(12)$ | $\hat{s}_1 \cdot \hat{s}_2$                                                         |
| $D(12)$      | $3(\hat{s}_1 \cdot \hat{r})(\hat{s}_2 \cdot \hat{r}) - (\hat{s}_1 \cdot \hat{s}_2)$ |

Los coeficientes radiales tienen una interpretación geométrica directa:
- $h_s(r)$: parte esférica (sin dependencia angular)
- $h_\Delta(r)$: correlación de alineación dipolo-dipolo
- $h_D(r)$: interacción tensorial a lo largo del eje intermolecular

### 2. Blum y Torruella (1972)

Introducen una expansión totalmente covariante usando invariantes rotacionales $\Phi^{mnl}_{\mu,\nu}$ construidos a partir de matrices de rotación de Wigner $D^m_{\mu'\mu}$ y símbolos $3j$ de Wigner:

$$\Phi^{mnl}_{\mu,\nu} = \sum_{\mu'\nu'\lambda'} \begin{pmatrix} m & n & l \\ \mu' & \nu' & \lambda' \end{pmatrix} D^m_{\mu',\mu}\, D^n_{\nu',\nu}\, D^l_{\lambda',0}$$

La expansión general es:

$$h(12) = \sum_{mnl\mu\nu} h^{mnl}_{\mu\nu}(r)\, \Phi^{mnl}_{\mu\nu}(12)$$

Para un sistema dipolar puro con partículas lineales, la suma se reduce a tres términos: $(000), (110), (112)$.

### 3. Patey (1977/1980)

Patey simplifica la notación de Blum introduciendo invariantes redefinidos con un prefactor explícito:

$$f^{mnl} = \frac{l!}{\begin{pmatrix} m & n & l \\ 0 & 0 & 0 \end{pmatrix}}$$

Esto permite definir $\Phi^{mnl}(12)$ de forma que coincida exactamente con las funciones de Wertheim:

$$\Phi^{110}(12) \equiv \Delta(12), \qquad \Phi^{112}(12) \equiv D(12)$$

La expansión queda:

$$h(12) = h^{000}(r) + h^{110}(r)\,\Phi^{110}(12) + h^{112}(r)\,\Phi^{112}(12)$$

---

## Conversiones entre representaciones

### Wertheim ↔ Patey

Gracias al prefactor $f^{mnl}$, los coeficientes son **numéricamente idénticos**:

| Patey        | Wertheim      | Descripción                                      |
| ------------ | ------------- | ------------------------------------------------ |
| $h^{000}(r)$ | $h_s(r)$      | Parte esférica / densidad                        |
| $h^{110}(r)$ | $h_\Delta(r)$ | Alineación dipolar ($\hat{s}_1 \cdot \hat{s}_2$) |
| $h^{112}(r)$ | $h_D(r)$      | Tensor dipolar ($D(12)$)                         |

> **Conclusión:** Si tienes datos en la convención de Wertheim, puedes usarlos directamente en el código sin ninguna conversión.

### Patey ↔ Blum-Torruella

Los coeficientes de Patey se obtienen de los de Blum multiplicando por el factor de escala:

$$\gamma^{mnl} = \frac{\left[(2m+1)(2n+1)\right]^{1/2}}{l!\, \begin{pmatrix} m & n & l \\ 0 & 0 & 0 \end{pmatrix}}$$

Concretamente, para los tres terminos dipolares:

| $(m,n,l)$ | $\begin{pmatrix} m & n & l \\ 0 & 0 & 0 \end{pmatrix}$ | $\gamma^{mnl}$         | Conversión                                                          |
| --------- | ------------------------------------------------------ | ---------------------- | ------------------------------------------------------------------- |
| $(0,0,0)$ | $1$                                                    | $1$                    | $h^{000}_\text{Patey} = h^{000}_\text{Blum}$                        |
| $(1,1,0)$ | $-1/\sqrt{3}$                                          | $3\sqrt{3}$            | $h^{110}_\text{Patey} = 3\sqrt{3}\, h^{110}_\text{Blum}$            |
| $(1,1,2)$ | $\sqrt{2/15}$                                          | $\frac{2\sqrt{30}}{2}$ | $h^{112}_\text{Patey} = \frac{2\sqrt{30}}{2}\, h^{112}_\text{Blum}$ |

> **Nota:** Los valores exactos de los símbolos 3j dependen de la convención de fase (Condon-Shortley). Verificar siempre con la referencia original de Blum y Torruella (1972).

---

## ¿Por qué Patey es la convención más útil?

1. **Compatibilidad con Wertheim**: Los coeficientes son idénticos, lo que facilita comparar con soluciones analíticas MSA (Wertheim 1971).
2. **Extensibilidad**: El esquema de índices $(m,n,l)$ permite agregar multipolos de orden superior (quadrupolos → $h^{224}$, $h^{123}$, etc.) sin cambiar el formalismo.
3. **Literatura dominante**: La mayoría de los papers de simulación y teoría de fluidos dipolares (Patey, Levesque, Weis, Gray & Gubbins) usan esta notación.

---

## La Representación Chi ($\chi$)

### Definición y origen

La representación $\chi$ es un conjunto de **combinaciones lineales de las proyecciones** en espacio de Fourier. Fue introducida por **L. Blum** como generalización de las funciones $h_+$ y $h_-$ originalmente definidas por Wertheim en su solución analítica MSA.

> **Ventaja principal:** Desacopla parcialmente la ecuación de Ornstein-Zernike en el espacio $k$, reduciendo la dimensionalidad del sistema de ecuaciones matriciales que hay que resolver. En lugar de invertir matrices de gran tamaño, el problema se descompone en subproblemas independientes más pequeños.

Por ejemplo, en cálculos RHNC para esferas duras dipolares con 29 proyecciones, la representación $\chi$ reduce la mayor matrix a invertir de $29 \times 29$ a solo $5 \times 5$.

### Fórmulas de conversión (Patey/Wertheim $\to$ $\chi$)

Las transformadas de Fourier (o Hankel) de cualquier función de correlación $f$ (puede ser $h$, $c$, o $\eta = h - c$) se combinan como:

$$f^\chi = \tilde{f}^{110} + \chi\, \tilde{f}^{112}$$

donde:

| Modo  | Fórmula en notación Patey              | Fórmula en notación Wertheim        |
| ----- | -------------------------------------- | ----------------------------------- |
| $f^0$ | $\tilde{f}^{110} + 2\,\tilde{f}^{112}$ | $\tilde{f}_\Delta + 2\,\tilde{f}_D$ |
| $f^1$ | $\tilde{f}^{110} - \tilde{f}^{112}$    | $\tilde{f}_\Delta - \tilde{f}_D$    |

La componente esférica $f^{000}$ (o $f_s$) **se desacopla completamente** de los términos anisotrópicos y se resuelve de forma independiente.

### Inversión: $\chi \to$ Patey/Wertheim

Dado que el sistema es:

$$f^0 = \tilde{f}^{110} + 2\tilde{f}^{112}, \qquad f^1 = \tilde{f}^{110} - \tilde{f}^{112}$$

se puede invertir directamente:

$$\tilde{f}^{110} = \frac{f^0 + 2f^1}{3}, \qquad \tilde{f}^{112} = \frac{f^0 - f^1}{3}$$

O equivalentemente en notación de Wertheim:

$$\tilde{f}_\Delta = \frac{f^0 + 2f^1}{3}, \qquad \tilde{f}_D = \frac{f^0 - f^1}{3}$$

### Diagrama de relaciones entre representaciones

```
                    ESPACIO REAL (r)
                    ────────────────
   Wertheim:    { hₛ(r),  h_Δ(r),  h_D(r) }
                        ↕ idénticos
   Patey:       { h⁰⁰⁰(r), h¹¹⁰(r), h¹¹²(r) }
                        ↕ × γ^{mnl}
   Blum-T.:     { h⁰⁰⁰_Blum, h¹¹⁰_Blum, h¹¹²_Blum }

                   ESPACIO FOURIER (k)
                   ───────────────────
   Patey:       { h̃⁰⁰⁰(k), h̃¹¹⁰(k), h̃¹¹²(k) }
                        ↕ combinación lineal
   Chi (Blum):  { h̃⁰⁰⁰(k),   h⁰(k),    h¹(k)  }
                  (esférico, desacoplado)    (χ=0, χ=1)
```

> **Nota para el código:** El solver actual trabaja con proyecciones en la convención de Patey en espacio real. Si se desea implementar la "Hat Transform" o el solver MSA analítico de Wertheim en espacio $k$, es conveniente hacer la transformación a la representación $\chi$ antes de la inversión matricial.

---

## Factor de Estructura y Ecuación de Ornstein-Zernike

### La ecuación de OZ en espacio de Fourier

Se define $\tilde{\eta}^{mnl}(k) = \tilde{h}^{mnl}(k) - \tilde{c}^{mnl}(k)$. La OZ para un sistema dipolar puro da **tres ecuaciones acopladas** (Patey 1977):

$$\tilde{\eta}^{000}(k) = \rho\, \tilde{c}^{000}(k)\, \tilde{h}^{000}(k)$$

$$\tilde{\eta}^{110}(k) = \frac{\rho}{3}\left[2\,\tilde{c}^{112}(k)\,\tilde{h}^{112}(k) + \tilde{c}^{110}(k)\,\tilde{h}^{110}(k)\right]$$

$$\tilde{\eta}^{112}(k) = \frac{\rho}{3}\left[\tilde{c}^{112}(k)\,\tilde{h}^{112}(k) + \tilde{c}^{112}(k)\,\tilde{h}^{110}(k) + \tilde{c}^{110}(k)\,\tilde{h}^{112}(k)\right]$$

La componente $h^{000}$ está **desacoplada**; $h^{110}$ y $h^{112}$ están **mutuamente acoplados**.

---

### Proyecciones del Factor de Estructura en notación Patey

El factor de estructura se define como:

$$S^{mnl}(k) = \delta_{mn0}\delta_{l0} + \rho\, \tilde{h}^{mnl}(k)$$

Las tres proyecciones en función de las correlaciones totales son directamente:

| Proyección | Fórmula |
|------------|---------|
| $S^{000}(k)$ | $1 + \rho\,\tilde{h}^{000}(k)$ |
| $S^{110}(k)$ | $\rho\,\tilde{h}^{110}(k)$ |
| $S^{112}(k)$ | $\rho\,\tilde{h}^{112}(k)$ |

En función sólo de las correlaciones directas $\tilde{c}^{mnl}$:

$$S^{000}(k) = \frac{1}{1 - \rho\,\tilde{c}^{000}(k)}$$

Para el bloque dipolar $\{110, 112\}$, se invierte el sistema acoplado. Definiendo el determinante:

$$\Delta(k) = \left(1 - \frac{\rho}{3}\tilde{c}^{110}\right)^2 - \left(\frac{\rho}{3}\tilde{c}^{112}\right)^2 - \frac{\rho}{3}\tilde{c}^{110}$$

las soluciones son:

$$S^{110}(k) = \frac{\rho\,\tilde{c}^{110}\left(1 - \frac{\rho}{3}\tilde{c}^{110}\right) + \frac{2\rho^2}{9}(\tilde{c}^{112})^2}{\Delta(k)}$$

$$S^{112}(k) = \frac{\rho\,\tilde{c}^{112}}{\Delta(k)}$$

---

### Factor de Estructura en la representación $\chi$

En la representación $\chi$ la OZ **se desacopla completamente**. Definiendo:

$$C^\chi(k) = \tilde{c}^{110}(k) + \chi\,\tilde{c}^{112}(k), \qquad \chi \in \{0,\, 1\}$$

cada modo satisface una OZ escalar independiente:

| Modo | Factor de estructura | Fórmula explícita |
|------|----------------------|-------------------|
| $\chi=0$ | $S^0(k)$ | $\dfrac{1}{1 - \dfrac{\rho}{3}(\tilde{c}^{110}+2\,\tilde{c}^{112})}$ |
| $\chi=1$ | $S^1(k)$ | $\dfrac{1}{1 + \dfrac{\rho}{3}(\tilde{c}^{110}-\tilde{c}^{112})}$ |

La componente esférica se trata independientemente: $S^{000}(k) = 1/(1 - \rho\,\tilde{c}^{000})$.

### Recuperando $S^{110}$ y $S^{112}$ desde la representación $\chi$

Con $S^0$ y $S^1$ calculados, se recuperan las proyecciones de Patey por inversión lineal:

$$S^{110}(k) = \frac{S^0(k) + 2\,S^1(k)}{3} - 1, \qquad S^{112}(k) = \frac{S^0(k) - S^1(k)}{3}$$

O equivalentemente en Wertheim:

$$1 + \rho\,\tilde{h}_\Delta = \frac{S^0 + 2S^1}{3}, \qquad \rho\,\tilde{h}_D = \frac{S^0 - S^1}{3}$$

> **Nota:** El modo $\chi=0$ generaliza $h_+$ de Wertheim y el modo $\chi=1$ generaliza $h_-$, cada uno satisfaciendo una OZ escalar con densidades efectivas $\rho/3$ y $-\rho/3$ respectivamente.

---

## Referencias

- Wertheim, M.S. (1971). *J. Chem. Phys.*, **55**, 4291.
- Blum, L. & Torruella, A.J. (1972). *J. Chem. Phys.*, **56**, 303.
- Patey, G.N. (1977). *Mol. Phys.*, **34**, 427.
- Patey, G.N., Levesque, D. & Weis, J.J. (1979). *Mol. Phys.*, **38**, 219.
- Gray, C.G. & Gubbins, K.E. (1984). *Theory of Molecular Fluids*, Vol. 1. Oxford.
