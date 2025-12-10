# OZE_c_solver - Solver de Ecuación de Ornstein-Zernike

Solver numérico para la ecuación de Ornstein-Zernike (OZ) aplicada a sistemas coloidales. Calcula el factor de estructura **S(k)** y la función de distribución radial **g(r)** usando las aproximaciones HNC (Hypernetted Chain) y RY (Rogers-Young).

## 🚀 Compilación Rápida

```bash
make           # Compilar el proyecto
make test      # Ejecutar prueba de ejemplo
```

## 📋 Requisitos

- **Compilador:** GCC (o compatible)
- **Bibliotecas:** GSL (GNU Scientific Library)
  ```bash
  # Ubuntu/Debian
  sudo apt-get install libgsl-dev
  
  # Fedora/RHEL
  sudo dnf install gsl-devel
  
  # macOS
  brew install gsl
  ```

## 🎯 Uso Básico

```bash
./build/facdes_solver --closure HNC --potential 13 \
                      --volfactor 0.3 --temp 1.0 \
                      --nodes 4096 --knodes 1024
```

### Parámetros Obligatorios:

| Parámetro     | Descripción                     | Ejemplo |
| ------------- | ------------------------------- | ------- |
| `--closure`   | Cierre termodinámico (HNC o RY) | `HNC`   |
| `--potential` | ID del potencial (1-13)         | `13`    |
| `--volfactor` | Fracción de volumen φ           | `0.3`   |
| `--temp`      | Temperatura T*                  | `1.0`   |
| `--nodes`     | Nodos espaciales                | `4096`  |
| `--knodes`    | Nodos en espacio k              | `1024`  |

### Parámetros Opcionales:

| Parámetro    | Descripción      | Default |
| ------------ | ---------------- | ------- |
| `--temp2`    | Temperatura T2*  | `1.0`   |
| `--lambda_a` | Lambda atractivo | `0.0`   |
| `--lambda_r` | Lambda repulsivo | `0.0`   |

## 🧪 Potenciales Disponibles

| ID  | Nombre                    | Ecuación               |
| --- | ------------------------- | ---------------------- |
| 1   | Inverse Power Law (IPL)   | U = T* (σ/r)^λ         |
| 2-3 | Truncated Lennard-Jones   | LJ truncado            |
| 4   | Double Yukawa             | Atractivo + Repulsivo  |
| 5   | Attractive Yukawa         | U ~ exp(-λr)/r         |
| 6   | Repulsive Yukawa          | U ~ exp(-λr)/r         |
| 7   | Hard Sphere (HS)          | U = ∞ (r<σ), 0 (r>σ)   |
| 8   | Shoulder Function         | Potencial tipo escalón |
| 9   | Down-Hill Function        | Lineal decreciente     |
| 10  | Gaussian Core Model       | U = T* exp(-(r/σ)²)    |
| 11  | Ramp (Step Function)      | U lineal tipo rampa    |
| 12  | Step Function (Soft Core) | U = E(1-r/σ)^n         |
| 13  | Hertzian Potential        | U = E(1-r/σ)^2.5       |

Ver ejemplos específicos: `./build/facdes_solver --potential <ID>`

## 📊 Archivos de Salida

El programa genera dos archivos principales:

- **`output/HNC_SdeK.dat`** (o `RY_SdeK.dat`): Factor de estructura S(k)
  ```
  k           S(k)
  0.000010    0.131638
  0.009823    0.131639
  ...
  ```

- **`output/HNC_GdeR.dat`** (o `RY_GdeR.dat`): Función de distribución radial g(r)
  ```
  r           g(r)
  0.000000    0.000000
  0.078125    0.000000
  1.015625    2.281424
  ...
  ```

## 📁 Estructura del Proyecto

```
OZE_c_solver/
├── src/                # Código fuente
│   ├── main.c
│   ├── facdes2Y.c
│   ├── math_aux.c
│   └── structures.c
├── include/            # Headers
│   ├── facdes2Y.h
│   ├── math_aux.h
│   └── structures.h
├── build/              # Ejecutable
│   └── facdes_solver
├── output/             # Archivos de salida .dat
├── examples/           # Scripts de ejemplo
├── docs/               # Documentación
├── Makefile
└── README.md
```

## 🛠️ Comandos Make Disponibles

```bash
make          # Compilar proyecto
make clean    # Limpiar archivos compilados
make cleanall # Limpiar todo (incluyendo .dat)
make test     # Ejecutar prueba de ejemplo
make help     # Mostrar ayuda
make install  # Instalar en /usr/local/bin
```

## 📖 Ejemplos

### Potencial Hertziano (n=2.5)
```bash
./build/facdes_solver --closure HNC --potential 13 \
                      --volfactor 0.3 --temp 1.0 \
                      --nodes 4096 --knodes 1024
```

### Hard Sphere
(Para esferas duras no es necesario el parámetro --temp, trabaja como una dump variable. Sin embargo, se debe proporcionar para cumplir con la sintaxis del programa)
```bash
./build/facdes_solver --closure HNC --potential 7 \
                      --volfactor 0.4 --temp 1.0 \
                      --nodes 4096 --knodes 1024
```

### Double Yukawa
```bash
./build/facdes_solver --closure HNC --potential 4 \
                      --volfactor 0.2 --temp 1.0 --temp2 0.5 \
                      --lambda_a 1.8 --lambda_r 5.0 \
                      --nodes 4096 --knodes 1024
```

## 📈 Visualización

```bash
# Con gnuplot
gnuplot
gnuplot> plot "output/HNC_GdeR.dat" with lines title "g(r)"
gnuplot> plot "output/HNC_SdeK.dat" with lines title "S(k)"
```

## 🔍 Teoría

El solver implementa el método de Ng para resolver iterativamente la ecuación de Ornstein-Zernike:

```
h(r) = c(r) + ρ ∫ c(|r-r'|) h(r') dr'
```

Usando relaciones de cierre:
- **HNC:** c(r) = exp(-βU(r) + γ(r)) - γ(r) - 1
- **RY:** Combinación de PY y HNC con parámetro α

## 📚 Referencias

- Hansen, J. P., & McDonald, I. R. (2013). *Theory of Simple Liquids*. Academic Press.
- Rogers, F. J., & Young, D. A. (1984). *Physical Review A*, 30(2), 999.

## 👤 Autor

Desarrollado por Ricardo Peredo Ortiz, Jonathan Josué Elisea Espinoza y el grupo de Materia fuera del equilibrio del Instituto de Física de la Universidad Autónoma de San Luis Potosí (UASLP). Se hace uso de herramientas de Inteligencia Artificial (Gemini 3) para el desarrollo de este proyecto. En particular, para la estructuración y organización del código.

## 📝 Licencia

MIT License
