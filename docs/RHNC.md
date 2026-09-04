### 1. Definición de Variables de Diferencia (Perturbación)
La teoría de referencia de Lado, adaptada por Fries y Patey, trabaja separando todas las funciones en una contribución del sistema de referencia esférico (subíndice \\(R\\), que para esferas duras dipolares es el sistema de esferas duras puras \\(HS\\)) y una perturbación anisotrópica (representada por el operador de diferencia \\(\Delta\\)):
\\[\Delta X(12) = X(12) - X_{R}(r)\\]

Dado que los canales anisotrópicos (como el \\(110\\), \\(112\\), etc.) tienen una contribución de referencia que es estrictamente cero (\\(c_{HS}^{mnl} = 0\\), \\(h_{HS}^{mnl} = 0\\), \\(\eta_{HS}^{mnl} = 0\\), etc.), las variables de diferencia se comportan de la siguiente manera en tu código:

*   **Para canales anisotrópicos (\\(mnl \neq 000\\)):**
    \\[\Delta c^{mnl}(r) = c^{mnl}(r)\\]
    \\[\Delta h^{mnl}(r) = h^{mnl}(r)\\]
    \\[\Delta \eta^{mnl}(r) = \eta^{mnl}(r)\\]
    \\[\Delta u^{mnl}(r) = u^{mnl}(r) = -\frac{\mu^2}{r^3} \quad (\text{solo para } 112 \text{ en dipolos puros})\\]

*   **Para el canal esférico (\\(000\\)):**
    Aquí es donde ocurre el único cambio real. Debes restar la solución exacta de las esferas duras puras (\\(c_{HS}(r)\\), \\(h_{HS}(r)\\), \\(\eta_{HS}(r)\\)) que usualmente alimentas como un *input* analítico desde la parametrización de Verlet-Weis:
    \\[\Delta c^{000}(r) = c^{000}(r) - c_{HS}(r)\\]
    \\[\Delta h^{000}(r) = h^{000}(r) - h_{HS}(r)\\]
    \\[\Delta \eta^{000}(r) = \eta^{000}(r) - \eta_{HS}(r)\\]
    \\[\Delta u^{000}(r) = u^{000}(r) - u_{HS}(r) = 0 \quad (\text{ya que ambos comparten el mismo núcleo duro})\\]

---

### 2. El tratamiento dentro del Núcleo Duro (\\(r < d\\))
Debido a la impenetrabilidad física, la función de distribución radial absoluta y la de referencia deben ser cero dentro del diámetro de contacto (\\(g(12) = g_{HS}(r) = 0\\)), lo que implica que \\(h(12) = h_{HS}(r) = -1\\). 
Sustituyendo esto en la ecuación de OZ, obtenemos una condición analítica exacta dentro del núcleo que es idéntica para todos los canales (incluyendo el esférico):
\\[\Delta c^{mnl}(r) = -\Delta \eta^{mnl}(r) \quad (\text{para } r < d)\\]

Esto simplifica enormemente el cómputo para \\(r < d\\), ya que no requiere evaluar exponenciales ni derivadas numéricas en la zona de solapamiento.

---

### 3. La Cerradura RHNC Fuera del Núcleo (\\(r > d\\)): El Gradiente
Para evitar la imposibilidad matemática de expandir directamente el término logarítmico no lineal de la cerradura, Fries y Patey derivan parcialmente la ecuación de RHNC respecto a la distancia radial \\(r\\). Al hacerlo, se obtiene una ecuación diferencial para la función de correlación directa de la perturbación:

\\[\frac{\partial \Delta c(12)}{\partial r} = -\Delta h(12) \frac{\partial \Delta W(12)}{\partial r} - h_{HS}(r)\frac{\partial \Delta W(12)}{\partial r} + \Delta h(12)\frac{\partial \ln g_{HS}(r)}{\partial r} - \beta \frac{\partial \Delta u(12)}{\partial r}\\]

Donde la función de potencial de fuerza media de la perturbación se define como:
\\[\Delta W(12) = -\Delta \eta(12) + \beta \Delta u(12) \implies \Delta W^{mnl}(r) = -\Delta \eta^{mnl}(r) + \beta \Delta u^{mnl}(r)\\]

Para implementar esto en tu código, debes notar que **la derivada total para cada canal \\((mnl)\\) se compone de 4 términos definidos en tu grilla**:

1.  **Término 1: Acoplamiento Anisotrópico (Anisotrópico \\(\times\\) Anisotrópico)**
    \\[T_1(12) = -\Delta h(12) \frac{\partial \Delta W(12)}{\partial r}\\]
    Este es el término verdaderamente no lineal. Como ya tienes QHNC o HNC, **debes usar exactamente el mismo motor de producto binario** que acopla los canales mediante símbolos 9-j de Wigner. Simplemente, en lugar de multiplicar las funciones absolutas, inyecta a tu función de multiplicación las proyecciones de \\(\Delta h^{mnl}(r)\\) y de las derivadas numéricas \\(\frac{\partial \Delta W^{mnl}(r)}{\partial r}\\):
    \\[T_1^{mnl}(r) = -\sum_{m_1 n_1 l_1} \sum_{m_2 n_2 l_2} P(m_1 n_1 l_1, m_2 n_2 l_2, m n l) \cdot \Delta h^{m_1 n_1 l_1}(r) \cdot \frac{\partial \Delta W^{m_2 n_2 l_2}(r)}{\partial r}\\]

2.  **Término 2: Acoplamiento de Referencia (Esférico \\(\times\\) Anisotrópico)**
    \\[T_2(12) = -h_{HS}(r)\frac{\partial \Delta W(12)}{\partial r}\\]
    Dado que la función de distribución de las esferas duras puras \\(h_{HS}(r)\\) es puramente esférica (solo tiene componente en el canal \\(000\\)), **este producto no acopla orientaciones angulares**. Para cada canal \\((mnl)\\), se calcula simplemente como un producto escalar radial directo (punto a punto en tu grilla, sin símbolos de Wigner):
    \\[T_2^{mnl}(r) = -h_{HS}(r) \frac{\partial \Delta W^{mnl}(r)}{\partial r}\\]

3.  **Término 3: Gradiente de la Referencia (Anisotrópico \\(\times\\) Esférico)**
    \\[T_3(12) = \Delta h(12)\frac{\partial \ln g_{HS}(r)}{\partial r}\\]
    La derivada de \\(\ln g_{HS}(r)\\) es también una función puramente esférica. Por ende, se evalúa directamente como otra multiplicación escalar radial para cada canal:
    \\[T_3^{mnl}(r) = \Delta h^{mnl}(r) \cdot \frac{\partial \ln g_{HS}(r)}{\partial r}\\]
    *(Nota: La derivada numérica de \\(\ln g_{HS}(r)\\) se precalcula una sola vez fuera del ciclo iterativo, ya que el sistema de referencia no cambia durante las iteraciones).*

4.  **Término 4: Derivada Analítica del Potencial**
    \\[T_4^{mnl}(r) = -\beta \frac{\partial \Delta u^{mnl}(r)}{\partial r}\\]
    Para dipolos puros, solo el canal \\(112\\) tiene potencial no nulo. Su derivada es analítica y exacta en toda la grilla fuera del núcleo:
    \\[T_4^{112}(r) = -\beta \frac{\partial}{\partial r}\left( -\frac{\mu^2}{r^3} \right) = -\frac{3\beta\mu^2}{r^4}\\]

---

### 4. Algoritmo de Integración de Regreso (Backwards Integration)
Una vez sumados los cuatro términos en cada punto radial para obtener la derivada total de cada canal, \\(\frac{\partial \Delta c^{mnl}(r)}{\partial r}\\):

1.  **Integración numérica:** Integra numéricamente hacia atrás utilizando diferencias finitas (trapecios o Simpson) partiendo desde el límite de tu caja de simulación \\(r = r_{max}\\) (donde asumes el límite asintótico \\(\Delta c^{mnl}(r_{max}) = 0\\)) hacia el contacto \\(r = d\\):
    \\[\Delta c^{mnl}(r)_{new} = \int_{r}^{r_{max}} \frac{\partial \Delta c^{mnl}(s)}{\partial s} ds\\]
2.  **Reconstrucción de las funciones absolutas:** Para poder alimentar el paso de la ecuación de Ornstein-Zernike, recuperas la función de correlación directa absoluta sumando la referencia en el canal esférico:
    *   **Canal esférico:** \\(c^{000}(r)_{new} = \Delta c^{000}(r)_{new} + c_{HS}(r)\\)
    *   **Canales anisotrópicos (\\(mnl \neq 000\\)):** \\(c^{mnl}(r)_{new} = \Delta c^{mnl}(r)_{new}\\)

¡Este nuevo conjunto de \\(c^{mnl}(r)_{new}\\) es el que mezclas (aplicando tus parámetros de mezcla \\(\alpha^{mnl}\\)) para comenzar la siguiente iteración!

---

### Resumen del ciclo iterativo de RHNC en tu código:
```
[ guesses iniciales c^{mnl}(r) ] 
       │
       ▼
 [ Calcular η^{mnl}(r) vía OZE en espacio k (Blum X-transforms) ]
       │
       ▼
 [ Restar referencia en canal 000 para obtener Δc, Δh, Δ\eta ]
       │
       ▼
 [ Aplicar cerradura: ]
 ├─► Para r < d: Δc^{mnl} = -Δ\eta^{mnl}
 └─► Para r > d: Evaluar los 4 términos de la derivada total,
                 integrar hacia atrás para hallar Δc^{mnl}_{new}
       │
       ▼
 [ Reconstruir c^{mnl}(r)_{new} sumando c_{HS}(r) al canal 000 ]
       │
       ▼
 [ Mezclar c^{(i+1)} = (1-α)c^{(i)} + α c_{new} y checar convergencia ]
```

¿Te gustaría que analicemos la fórmula exacta o el código de Verlet-Weis para inicializar correctamente las funciones de referencia \\(g_{HS}(r)\\) y \\(c_{HS}(r)\\) en tu malla radial, o prefieres revisar cómo acoplar las derivadas de la perturbación \\(\Delta W^{mnl}\\) usando tu motor de símbolos 9-j?
