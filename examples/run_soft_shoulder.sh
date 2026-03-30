#!/bin/bash
# Example script: Soft Shoulder Potential (ID 16)

# Change to the build directory (assuming the script is in examples/)
cd "$(dirname "$0")/.." || exit 1

# Check if the executable exists, try to build if not
if [ ! -f build/facdes_solver ]; then
    echo "Executable build/facdes_solver not found. Running make..."
    make
fi

# Run the solver for Soft Shoulder Potential
# parameters: phi=0.2, T=1.0, lambda=2.0, alpha=5.0
./build/facdes_solver --closure HNC --potential 16 \
                --volfactor 0.2 --temp 1.0 \
                --lambda_a 2.0 --lambda_r 5.0 \
                --nodes 2048 --knodes 1024

echo ""
echo "========================================="
echo " Files generated in output/:"
echo "========================================="
ls -lh output/*.dat 2>/dev/null || echo "No .dat files generated"
echo ""
echo "To visualize with gnuplot:"
echo "  gnuplot> plot \"output/HNC_GdeR.dat\" with lines title \"g(r)\""
echo "  gnuplot> plot \"output/HNC_SdeK.dat\" with lines title \"S(k)\""
