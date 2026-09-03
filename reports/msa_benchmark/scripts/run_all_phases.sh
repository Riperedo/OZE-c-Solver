#!/usr/bin/env bash
# ==============================================================================
# Master Execution Script: Three-Phase MSA Benchmark Suite
# ==============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

cd "${ROOT_DIR}"

echo "========================================================================"
echo "1. Building C Solver..."
echo "========================================================================"
make -j4

echo ""
echo "========================================================================"
echo "2. Running Phase 1: Standard Cold-Start Solver (No Initializer)..."
echo "========================================================================"
python3 reports/msa_benchmark/scripts/run_phase1_cold_start.py

echo ""
echo "========================================================================"
echo "3. Running Phase 2: Temperature Continuation Ramps (Annealing)..."
echo "========================================================================"
python3 reports/msa_benchmark/scripts/run_phase2_ramps.py

echo ""
echo "========================================================================"
echo "4. Running Phase 3: Direct Analytical S(k) Warm-Start..."
echo "========================================================================"
python3 reports/msa_benchmark/scripts/run_phase3_ana_init.py

echo ""
echo "========================================================================"
echo "5. Generating Publication Figures (Figures 1 to 7)..."
echo "========================================================================"
gnuplot reports/msa_benchmark/scripts/generate_plots.gp
gnuplot reports/msa_benchmark/scripts/generate_fig5.gp
python3 reports/msa_benchmark/scripts/generate_heatmap.py
python3 reports/msa_benchmark/scripts/generate_fig7_ana_init.py

echo ""
echo "========================================================================"
echo "6. Compiling LaTeX Academic Report..."
echo "========================================================================"
cd reports/msa_benchmark
pdflatex -interaction=nonstopmode msa_benchmark_report.tex > /dev/null
pdflatex -interaction=nonstopmode msa_benchmark_report.tex > /dev/null
cd "${ROOT_DIR}"

echo ""
echo "========================================================================"
echo "✓ All benchmarks completed and report compiled successfully!"
echo "Report PDF: reports/msa_benchmark/msa_benchmark_report.pdf"
echo "========================================================================"
