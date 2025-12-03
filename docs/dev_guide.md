# Guía de Desarrollador - OZE_c_solver

Esta guía está destinada a desarrolladores que deseen modificar, extender o entender la estructura interna del código.

## 1. Estructura del Proyecto

```
OZE_c_solver/
├── src/                # Código fuente (.c)
│   ├── main.c          # Punto de entrada, parseo de argumentos
│   ├── facdes2Y.c      # Interfaz de alto nivel, gestión de memoria
│   ├── structures.c    # Núcleo del solver (potenciales, Ng, cierres)
│   └── math_aux.c      # Funciones matemáticas (FFT, integrales)
├── include/            # Archivos de cabecera (.h)
├── build/              # Archivos objeto y ejecutable
├── output/             # Archivos de salida generados
├── docs/               # Documentación
└── Makefile            # Sistema de compilación
```

## 2. Flujo de Ejecución

1.  **`main.c`**: Lee los argumentos de la línea de comandos, valida la entrada y llama a las funciones de `facdes2Y.c` (e.g., `sk_HNC`, `gr_RY`).
2.  **`facdes2Y.c`**:
    - Reserva memoria para los vectores de trabajo.
    - Llama a `facdes2YFunc` (función wrapper).
    - Llama a `input` (en `structures.c`) para inicializar potenciales.
    - Llama a `OZ2` (en `structures.c`) para resolver la ecuación.
    - Interpola los resultados y los escribe en disco.
3.  **`structures.c`**:
    - `input()`: Inicializa la malla y llama a `POT()`.
    - `POT()`: Calcula el potencial $U(r)$ y su derivada.
    - `OZ2()`: Implementa el bucle principal del método de Ng y la rampa de densidad.
    - `Ng()`: Algoritmo de aceleración de convergencia.
    - `closrel()`: Aplica la relación de cierre (HNC, RY, PY).

## 3. Cómo Añadir un Nuevo Potencial

Para añadir un nuevo potencial de interacción (digamos, ID `14`):

1.  Abra `src/structures.c`.
2.  Busque la función `POT`.
3.  Añada un nuevo `case 14:` dentro del `switch(potentialID)`.
4.  Implemente el cálculo de `U[i*ncols + k]` (potencial) y `Up[i*ncols + k]` (derivada $-dU/dr \cdot r$ o similar, verifique consistencia con otros casos).
    - **Nota**: `Up` se usa para el cálculo de la presión virial.
5.  Añada la descripción en `PotentialName` (al final de `structures.c`).
6.  (Opcional) Actualice `display_potential_options` en `src/main.c` para que aparezca en la ayuda.

### Ejemplo de Estructura
```c
case 14:
    // Inicializar parámetros (E, z, etc.) usando especie1.temperature, etc.
    // ...
    for (k = 0; k < ncols; k++) {
        for (i = 0; i < nrows; i++) {
            // Calcular r[i]
            // U[i*ncols + k] = ...
            // Up[i*ncols + k] = ...
        }
    }
    break;
```

## 4. Cómo Añadir un Nuevo Cierre

Para añadir una nueva relación de cierre (e.g., Martynov-Sarkisov):

1.  Abra `src/structures.c`.
2.  Busque la función `closrel`.
3.  Añada un nuevo `case` en el `switch(closureID)`.
4.  Implemente la relación $c(r) = f(h(r), U(r))$.
5.  Actualice `main.c` para aceptar el nuevo string en el argumento `--closure`.

## 5. Estilo de Código

- **Indentación**: 4 espacios.
- **Nombres**: CamelCase para funciones (`facdes2YFunc`), minúsculas para variables locales.
- **Memoria**: Siempre liberar (`free`) la memoria reservada con `malloc`.
