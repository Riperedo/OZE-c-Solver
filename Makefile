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
SOURCES = $(SRC_DIR)/main.c $(SRC_DIR)/facdes2Y.c $(SRC_DIR)/math_aux.c $(SRC_DIR)/structures.c
HEADERS = $(INC_DIR)/facdes2Y.h $(INC_DIR)/math_aux.h $(INC_DIR)/structures.h
OBJECTS = $(BUILD_DIR)/main.o $(BUILD_DIR)/facdes2Y.o $(BUILD_DIR)/math_aux.o $(BUILD_DIR)/structures.o
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

# Instalar (copiar ejecutable a /usr/local/bin)
install: $(TARGET)
	@echo "Instalando facdes_solver..."
	sudo cp $(TARGET) /usr/local/bin/
	@echo "$(GREEN)✓ Instalado en /usr/local/bin/facdes_solver$(NC)"

# Desinstalar
uninstall:
	@echo "Desinstalando facdes_solver..."
	sudo rm -f /usr/local/bin/facdes_solver
	@echo "$(GREEN)✓ Desinstalado!$(NC)"

.PHONY: all clean cleanall dirs test help install uninstall
