# Makefile for OZE_c_solver
# Ornstein-Zernike Equation Solver (Spherical & Non-Spherical Closures)

# Compiler and Flags
CC = gcc
CFLAGS = -Wall -O2 -Iinclude
LIBS = -lgsl -lgslcblas -lm

# Directories
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build
OUT_DIR = output

# Source and Object Files
SOURCES = $(SRC_DIR)/main.c $(SRC_DIR)/facdes2Y.c $(SRC_DIR)/math_aux.c $(SRC_DIR)/structures.c $(SRC_DIR)/structures_nonspherical.c $(SRC_DIR)/closures_nonspherical.c $(SRC_DIR)/solver_dipolar.c $(SRC_DIR)/solver_mode2.c
HEADERS = $(INC_DIR)/facdes2Y.h $(INC_DIR)/math_aux.h $(INC_DIR)/structures.h $(INC_DIR)/structures_nonspherical.h
OBJECTS = $(BUILD_DIR)/main.o $(BUILD_DIR)/facdes2Y.o $(BUILD_DIR)/math_aux.o $(BUILD_DIR)/structures.o $(BUILD_DIR)/structures_nonspherical.o $(BUILD_DIR)/closures_nonspherical.o $(BUILD_DIR)/solver_dipolar.o $(BUILD_DIR)/solver_mode2.o
TARGET = $(BUILD_DIR)/facdes_solver

# Terminal Colors
GREEN = \033[0;32m
NC = \033[0m # No Color

# Default Rule
all: $(TARGET)
	@echo "$(GREEN)✓ Build successful!$(NC)"
	@echo "Executable: $(TARGET)"

# Executable Linking Rule
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(CC) $(OBJECTS) $(LIBS) -o $(TARGET)

# Object Files Compilation Rule
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	@mkdir -p $(BUILD_DIR)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# Clean Compiled Artifacts
clean:
	@echo "Cleaning compiled objects..."
	rm -f $(BUILD_DIR)/*.o $(TARGET)
	@echo "$(GREEN)✓ Clean complete!$(NC)"

# Clean All (including output .dat files)
cleanall: clean
	@echo "Cleaning output data files..."
	rm -f $(OUT_DIR)/*.dat
	@echo "$(GREEN)✓ Clean complete (including .dat outputs)!$(NC)"

# Create Directory Structure
dirs:
	@mkdir -p $(SRC_DIR) $(INC_DIR) $(BUILD_DIR) $(OUT_DIR) examples docs reports

# Fast Spherical Test (Hertzian)
test: $(TARGET)
	@echo "Running spherical test with Hertzian potential..."
	@./$(TARGET) --closure HNC --potential 13 --volfactor 0.3 --temp 1.0 --nodes 2048 --knodes 512
	@echo "$(GREEN)✓ Spherical test completed successfully!$(NC)"
	@echo "Output generated in $(OUT_DIR)/"

# Fast Dipolar Closures Test (MSA, LHNC, QHNC, RHNC)
test-dipolar: $(TARGET)
	@echo "Running dipolar closures test suite (MSA, LHNC, QHNC, RHNC)..."
	@echo "  -> Testing MSA..."
	@./$(TARGET) --closure MSA --potential 14 --volfactor 0.418879 --temp 1.0 --dipole 1.0 --nodes 2048 --rmax 15.0 > /dev/null
	@echo "  -> Testing LHNC..."
	@./$(TARGET) --closure LHNC --potential 14 --volfactor 0.418879 --temp 1.0 --dipole 1.0 --nodes 2048 --rmax 15.0 > /dev/null
	@echo "  -> Testing QHNC..."
	@./$(TARGET) --closure QHNC --potential 14 --volfactor 0.418879 --temp 1.0 --dipole 1.0 --nodes 2048 --rmax 15.0 > /dev/null
	@echo "  -> Testing RHNC..."
	@./$(TARGET) --closure RHNC --potential 14 --volfactor 0.418879 --temp 1.0 --dipole 1.0 --nodes 2048 --rmax 15.0 > /dev/null
	@echo "$(GREEN)✓ All dipolar closures (MSA, LHNC, QHNC, RHNC) verified successfully!$(NC)"

# Complete Test Suite
test-all: test test-dipolar
	@echo "$(GREEN)✓ Full test suite passed successfully!$(NC)"

# Help Menu
help:
	@echo "Makefile for OZE_c_solver - Ornstein-Zernike Equation Solver"
	@echo ""
	@echo "Usage:"
	@echo "  make               - Build the executable"
	@echo "  make clean         - Remove object files and executable"
	@echo "  make cleanall      - Remove objects and generated .dat outputs"
	@echo "  make test          - Run fast spherical HNC test"
	@echo "  make test-dipolar  - Run fast dipolar test suite (MSA, LHNC, QHNC, RHNC)"
	@echo "  make test-all      - Run all verification tests"
	@echo "  make mc-all        - Run full Monte Carlo benchmark pipeline (benchmarks, plots, report)"
	@echo "  make report-mc     - Compile the Monte Carlo comparison academic PDF report"
	@echo "  make benchmark-all - Run Wertheim MSA analytical benchmark pipeline"
	@echo "  make report        - Compile the MSA benchmark academic PDF report"
	@echo "  make help          - Display this help message"

# MSA Benchmark and Report Targets
benchmark-phase1: $(TARGET)
	@echo "Running Phase 1: Standard Cold-Start Benchmark..."
	python3 reports/msa_benchmark/scripts/run_phase1_cold_start.py

benchmark-phase2: $(TARGET)
	@echo "Running Phase 2: Temperature Continuation Ramps..."
	python3 reports/msa_benchmark/scripts/run_phase2_ramps.py

benchmark-phase3: $(TARGET)
	@echo "Running Phase 3: Direct Analytical S(k) Warm-Start..."
	python3 reports/msa_benchmark/scripts/run_phase3_ana_init.py

benchmark-plots:
	@echo "Generating Publication Figures (Figures 1 to 7)..."
	gnuplot reports/msa_benchmark/scripts/generate_plots.gp
	gnuplot reports/msa_benchmark/scripts/generate_fig5.gp
	python3 reports/msa_benchmark/scripts/generate_heatmap.py
	python3 reports/msa_benchmark/scripts/generate_fig7_ana_init.py

report:
	@echo "Compiling LaTeX Academic Report (MSA Benchmark)..."
	@cd reports/msa_benchmark && pdflatex -interaction=nonstopmode msa_benchmark_report.tex > /dev/null && pdflatex -interaction=nonstopmode msa_benchmark_report.tex > /dev/null
	@echo "$(GREEN)✓ Report compiled: reports/msa_benchmark/msa_benchmark_report.pdf$(NC)"

benchmark-all:
	@bash reports/msa_benchmark/scripts/run_all_phases.sh

# Monte Carlo Simulation Benchmark Targets (Fries & Patey 1985 + Patey, Levesque & Weis 1979)
benchmark-mc: $(TARGET)
	@echo "Running Fries & Patey (1985) Benchmark Suite (MSA, LHNC, QHNC, RHNC)..."
	python3 reports/monte_carlo_benchmark/scripts/run_mc_benchmarks.py

plots-mc:
	@echo "Generating Fries & Patey (1985) Publication Figures (Figures 1 to 6)..."
	python3 reports/monte_carlo_benchmark/scripts/generate_mc_plots.py

benchmark-plw: $(TARGET)
	@echo "Running Patey, Levesque & Weis (1979) Benchmark Suite (MSA, LHNC, QHNC, RHNC)..."
	python3 reports/monte_carlo_benchmark/scripts/run_plw_benchmarks.py

plots-plw:
	@echo "Generating Patey, Levesque & Weis (1979) Publication Figures..."
	python3 reports/monte_carlo_benchmark/scripts/generate_plw_plots.py

report-mc: plots-mc plots-plw
	@echo "Compiling Comprehensive Monte Carlo Academic Report..."
	@cd reports/monte_carlo_benchmark && pdflatex -interaction=nonstopmode monte_carlo_benchmark_report.tex > /dev/null && pdflatex -interaction=nonstopmode monte_carlo_benchmark_report.tex > /dev/null
	@echo "$(GREEN)✓ Monte Carlo Report compiled: reports/monte_carlo_benchmark/monte_carlo_benchmark_report.pdf$(NC)"

# Percus-Yevick Hard Sphere Benchmark Targets (Wertheim-Thiele Exact Analytical Solution)
benchmark-py: $(TARGET)
	@echo "Running Percus-Yevick Hard Sphere Benchmark Suite across 7 packing fractions..."
	python3 reports/py_hard_sphere_benchmark/scripts/run_py_benchmark.py

plots-py:
	@echo "Generating Percus-Yevick Publication Figures (Figures 1 to 5)..."
	python3 reports/py_hard_sphere_benchmark/scripts/generate_py_plots.py

report-py: plots-py
	@echo "Compiling Percus-Yevick Hard Sphere Academic Report..."
	@cd reports/py_hard_sphere_benchmark && pdflatex -interaction=nonstopmode py_hard_sphere_report.tex > /dev/null && pdflatex -interaction=nonstopmode py_hard_sphere_report.tex > /dev/null
	@echo "$(GREEN)✓ Percus-Yevick Report compiled: reports/py_hard_sphere_benchmark/py_hard_sphere_report.pdf$(NC)"

py-all: benchmark-py plots-py report-py

.PHONY: all clean cleanall dirs test test-dipolar test-all help benchmark-phase1 benchmark-phase2 benchmark-phase3 benchmark-plots report benchmark-all benchmark-mc plots-mc benchmark-plw plots-plw report-mc mc-all benchmark-py plots-py report-py py-all


