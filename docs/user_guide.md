# Guía de Usuario - OZE_c_solver

Bienvenido a la guía de usuario del solver de la ecuación de Ornstein-Zernike. Este documento explica cómo compilar, ejecutar y configurar el programa para diferentes sistemas coloidales.

## 1. Instalación y Compilación

### Requisitos Previos
Asegúrese de tener instalado `gcc` y la librería `GSL` (GNU Scientific Library).

```bash
# Ubuntu/Debian
sudo apt-get install libgsl-dev build-essential

# Fedora
sudo dnf install gsl-devel
```

### Compilación
Para compilar el proyecto, simplemente ejecute `make` en la raíz del directorio:

```bash
make
```

Esto generará el ejecutable `build/facdes_solver`.

## 2. Ejecución Básica

El programa se ejecuta desde la línea de comandos. La sintaxis general es:

```bash
./build/facdes_solver --closure [TIPO] --potential [ID] --volfactor [PHI] --temp [T] --nodes [N] --knodes [Nk]
```

### Argumentos Obligatorios

| Argumento     | Descripción                                                                             | Ejemplo |
| :------------ | :-------------------------------------------------------------------------------------- | :------ |
| `--closure`   | Relación de cierre: `HNC` (Hypernetted Chain) o `RY` (Rogers-Young).                    | `HNC`   |
| `--potential` | ID numérico del potencial de interacción (ver sección 3).                               | `13`    |
| `--volfactor` | Fracción de volumen del sistema ($\phi$).                                               | `0.3`   |
| `--temp`      | Temperatura reducida ($T^*$) o parámetro de energía.                                    | `1.0`   |
| `--nodes`     | Número de puntos en la malla espacial (espacio real $r$). Se recomienda potencias de 2. | `4096`  |
| `--knodes`    | Número de puntos en el espacio recíproco ($k$) para el archivo de salida.               | `1024`  |

### Argumentos Opcionales

| Argumento    | Descripción                                                        | Default |
| :----------- | :----------------------------------------------------------------- | :------ |
| `--temp2`    | Segunda temperatura o parámetro de ancho para ciertos potenciales. | `1.0`   |
| `--lambda_a` | Parámetro de alcance atractivo o exponente.                        | `0.0`   |
| `--lambda_r` | Parámetro de alcance repulsivo.                                    | `0.0`   |

## 3. Catálogo de Potenciales

A continuación se detallan los potenciales disponibles y sus parámetros específicos.

### 1. Inverse Power Law (IPL)
$$ U(r) = \epsilon (\sigma/r)^\lambda $$
- **ID**: `1`
- **Parámetros**:
    - `--temp`: $\epsilon$ (Energía)
    - `--lambda_a`: $\lambda$ (Exponente)

### 2. Lennard-Jones Truncado (Repulsivo)
Solo la parte repulsiva del potencial LJ (WCA).
- **ID**: `2`
- **Parámetros**:
    - `--temp`: Temperatura reducida $T^*$

### 4. Double Yukawa
Combinación de una parte atractiva y una repulsiva.
$$ U(r) = -K_a \frac{e^{-\lambda_a r}}{r} + K_r \frac{e^{-\lambda_r r}}{r} $$
- **ID**: `4`
- **Parámetros**:
    - `--temp`: Inverso de la intensidad atractiva ($1/K_a$)
    - `--temp2`: Inverso de la intensidad repulsiva ($1/K_r$)
    - `--lambda_a`: Alcance atractivo $\lambda_a$
    - `--lambda_r`: Alcance repulsivo $\lambda_r$

### 7. Hard Sphere (Esferas Duras)
- **ID**: `7`
- **Parámetros**:
    - `--temp`: Valor dummy (no afecta física, usar 1.0)

### 13. Hertzian Potential
Modelo de esferas blandas elásticas.
$$ U(r) = \epsilon (1 - r/\sigma)^{5/2} \quad \text{si } r < \sigma $$
- **ID**: `13`
- **Parámetros**:
    - `--temp`: Energía $\epsilon$

*(Para ver la lista completa, ejecute `./build/facdes_solver` sin argumentos)*

## 4. Archivos de Salida

Los resultados se guardan en la carpeta `output/`.

- **`HNC_SdeK.dat`** / **`RY_SdeK.dat`**: Factor de estructura estático $S(k)$.
    - Columna 1: Vector de onda $k$
    - Columna 2: $S(k)$
- **`HNC_GdeR.dat`** / **`RY_GdeR.dat`**: Función de distribución radial $g(r)$.
    - Columna 1: Distancia $r$
    - Columna 2: $g(r)$

## 5. Ejemplos Prácticos

### Ejemplo 1: Esferas Duras (Hard Spheres)
Simulación de un líquido de esferas duras a densidad media.

```bash
./build/facdes_solver --closure PY --potential 7 --volfactor 0.4 --temp 1.0 --nodes 2048 --knodes 512
```
*(Nota: El código usa HNC o RY, para HS puro RY suele ser mejor o equivalente a PY si alpha se ajusta)*

### Ejemplo 2: Potencial Hertziano
Simulación de coloides blandos.

```bash
./build/facdes_solver --closure HNC --potential 13 --volfactor 0.3 --temp 1.0 --nodes 4096 --knodes 1024
```

### Visualización Rápida con Gnuplot
```bash
gnuplot -e "plot 'output/HNC_SdeK.dat' w l title 'S(k)'; pause -1"
```
