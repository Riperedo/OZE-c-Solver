#!/bin/bash
# Script de ejemplo: Potencial Yukawa

# Cambiar al directorio build
cd "$(dirname "$0")/../build" || exit 1

# Ejecutar el solver
./facdes_solver --closure HNC --potential 13 \
                --volfactor 0.8 --temp 0.5 \
                --nodes 2048 --knodes 512

echo ""
echo "========================================="
echo " Archivos generados:"
echo "========================================="
ls -lh ../output/*.dat 2>/dev/null || echo "No se generaron archivos .dat"
echo ""
echo "Para visualizar con gnuplot:"
echo "  gnuplot"
echo "  gnuplot> plot \"../output/HNC_GdeR.dat\" with lines title \"g(r)\""
echo "  gnuplot> plot \"../output/HNC_SdeK.dat\" with lines title \"S(k)\""
