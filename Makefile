# Makefile para API_HNC
# Solver de Ecuación de Ornstein-Zernike

# Compilador y flags
CC = gcc
CFLAGS = -Wall -O2 -Iinclude
LIBS = -lgsl -lgslcblas -lm

# Directorios
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build
OUT_DIR = output

# Archivos fuente y objeto
SOURCES = $(SRC_DIR)/main.c $(SRC_DIR)/facdes2Y.c $(SRC_DIR)/math_aux.c $(SRC_DIR)/structures.c $(SRC_DIR)/structures_nonspherical.c $(SRC_DIR)/closures_nonspherical.c $(SRC_DIR)/solver_dipolar.c $(SRC_DIR)/solver_mode2.c
HEADERS = $(INC_DIR)/facdes2Y.h $(INC_DIR)/math_aux.h $(INC_DIR)/structures.h $(INC_DIR)/structures_nonspherical.h
OBJECTS = $(BUILD_DIR)/main.o $(BUILD_DIR)/facdes2Y.o $(BUILD_DIR)/math_aux.o $(BUILD_DIR)/structures.o $(BUILD_DIR)/structures_nonspherical.o $(BUILD_DIR)/closures_nonspherical.o $(BUILD_DIR)/solver_dipolar.o $(BUILD_DIR)/solver_mode2.o
TARGET = $(BUILD_DIR)/facdes_solver

# Colores para output
GREEN = \033[0;32m
NC = \033[0m # No Color

# Regla por defecto
all: $(TARGET)
	@echo "$(GREEN)✓ Compilación exitosa!$(NC)"
	@echo "Ejecutable: $(TARGET)"

# Regla para el ejecutable
$(TARGET): $(OBJECTS)
	@echo "Enlazando $(TARGET)..."
	$(CC) $(OBJECTS) $(LIBS) -o $(TARGET)

# Reglas para archivos objeto
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	@mkdir -p $(BUILD_DIR)
	@echo "Compilando $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# Limpiar archivos compilados
clean:
	@echo "Limpiando archivos compilados..."
	rm -f $(BUILD_DIR)/*.o $(TARGET)
	@echo "$(GREEN)✓ Limpieza completa!$(NC)"

# Limpiar todo (incluyendo salidas)
cleanall: clean
	@echo "Limpiando archivos de salida..."
	rm -f $(OUT_DIR)/*.dat
	@echo "$(GREEN)✓ Limpieza completa (incluyendo .dat)!$(NC)"

# Crear directorios necesarios
dirs:
	@mkdir -p $(SRC_DIR) $(INC_DIR) $(BUILD_DIR) $(OUT_DIR) examples docs

# Ejecutar ejemplo con Hertzian
test: $(TARGET)
	@echo "Ejecutando prueba con potencial Hertziano..."
	@./$(TARGET) --closure HNC --potential 13 --volfactor 0.3 --temp 1.0 --nodes 2048 --knodes 512
	@echo "$(GREEN)✓ Prueba completada!$(NC)"
	@echo "Archivos generados en $(OUT_DIR)/"

# Mostrar ayuda
help:
	@echo "Makefile para API_HNC - Solver de Ecuación de Ornstein-Zernike"
	@echo ""
	@echo "Uso:"
	@echo "  make          - Compilar el proyecto"
	@echo "  make clean    - Limpiar archivos compilados"
	@echo "  make cleanall - Limpiar todo (incluyendo .dat)"
	@echo "  make test     - Ejecutar prueba de ejemplo"
	@echo "  make help     - Mostrar esta ayuda"
	@echo ""
	@echo "Ejemplo de ejecución manual:"
	@echo "  $(TARGET) --closure HNC --potential 13 --volfactor 0.3 --temp 1.0 --nodes 2048 --knodes 512"

# Benchmark and Report Targets
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

# Monte Carlo Simulation Benchmark Targets (Fries & Patey JCP 1985)
benchmark-mc: $(TARGET)
	@echo "Running Monte Carlo Benchmark Suite (MSA, LHNC, QHNC)..."
	python3 reports/monte_carlo_benchmark/scripts/run_mc_benchmarks.py

plots-mc:
	@echo "Generating Monte Carlo Publication Figures (Figures 1 to 6)..."
	python3 reports/monte_carlo_benchmark/scripts/generate_mc_plots.py

report-mc: plots-mc
	@echo "Compiling Monte Carlo Academic Report..."
	@cd reports/monte_carlo_benchmark && pdflatex -interaction=nonstopmode monte_carlo_benchmark_report.tex > /dev/null && pdflatex -interaction=nonstopmode monte_carlo_benchmark_report.tex > /dev/null
	@echo "$(GREEN)✓ Monte Carlo Report compiled: reports/monte_carlo_benchmark/monte_carlo_benchmark_report.pdf$(NC)"

mc-all: benchmark-mc plots-mc report-mc

.PHONY: all clean cleanall dirs test help install uninstall benchmark-phase1 benchmark-phase2 benchmark-phase3 benchmark-plots report benchmark-all benchmark-mc plots-mc report-mc mc-all
