#!/bin/bash
# Script para ejecutar clean.sh y run.sh en todas las carpetas pump=*

# Colores para mejor visualización
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Ejecutando simulaciones en carpetas pump${NC}"
echo -e "${BLUE}========================================${NC}\n"

# Directorio base (donde está este script)
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Contador de carpetas procesadas
count=0

# Recorrer todas las carpetas que empiezan con "pump="
for dir in "$BASE_DIR"/pump=*/; do
    # Verificar que el directorio existe
    if [ ! -d "$dir" ]; then
        continue
    fi
    
    # Obtener el nombre de la carpeta
    folder_name=$(basename "$dir")
    
    echo -e "${YELLOW}================================================${NC}"
    echo -e "${YELLOW}Procesando: $folder_name${NC}"
    echo -e "${YELLOW}================================================${NC}"
    
    # Entrar a la carpeta
    cd "$dir" || {
        echo -e "${RED}Error: No se pudo entrar a la carpeta $folder_name${NC}"
        continue
    }
    
    # Ejecutar clean.sh
    if [ -f "clean.sh" ]; then
        echo -e "${GREEN}→ Ejecutando ./clean.sh${NC}"
        chmod +x clean.sh
        ./clean.sh
    else
        echo -e "${RED}⚠️  Advertencia: clean.sh no encontrado en $folder_name${NC}"
    fi
    
    # Ejecutar run.sh con nohup en background
    if [ -f "run.sh" ]; then
        echo -e "${GREEN}→ Ejecutando nohup ./run.sh &${NC}"
        chmod +x run.sh
        nohup ./run.sh > nohup.out 2>&1 &
        pid=$!
        echo -e "${GREEN}✓ Proceso iniciado con PID: $pid${NC}"
    else
        echo -e "${RED}⚠️  Advertencia: run.sh no encontrado en $folder_name${NC}"
    fi
    
    echo ""
    
    # Volver al directorio base
    cd "$BASE_DIR"
    
    ((count++))
done

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Resumen:${NC}"
echo -e "${BLUE}  - Carpetas procesadas: $count${NC}"
echo -e "${BLUE}  - Procesos en ejecución en background${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "\n${GREEN}Para ver los procesos en ejecución:${NC}"
echo -e "  ps aux | grep run.sh"
echo -e "\n${GREEN}Para ver el output de una simulación:${NC}"
echo -e "  tail -f pump=X.X/nohup.out"
echo -e "\n${GREEN}Para detener todos los procesos:${NC}"
echo -e "  pkill -f run.sh"
