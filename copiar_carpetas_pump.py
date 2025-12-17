#!/usr/bin/env python3
"""
Script para copiar la carpeta "Main Final" con diferentes valores de pump.
Crea carpetas pump=0.1, pump=0.2, ..., pump=1.0 y modifica el archivo input/pump en cada una.
"""

import shutil
import os
from pathlib import Path


def main():
    # Directorio base (donde está este script)
    base_dir = Path(__file__).parent
    
    # Carpeta fuente
    source_folder = base_dir / "Main Final"
    
    # Verificar que la carpeta fuente existe
    if not source_folder.exists():
        print(f"Error: La carpeta '{source_folder}' no existe.")
        return
    
    # Valores de pump: 0.1, 0.2, ..., 1.0
    pump_values = [round(i * 0.01, 2) for i in range(1, 11)]
    
    print(f"Copiando carpeta '{source_folder.name}' con valores de pump: {pump_values}\n")
    
    for pump_value in pump_values:
        # Nombre de la nueva carpeta
        new_folder_name = f"pump={pump_value}"
        new_folder_path = base_dir / new_folder_name
        
        # Copiar la carpeta
        print(f"Creando carpeta: {new_folder_name}")
        if new_folder_path.exists():
            print(f"  ⚠️  La carpeta '{new_folder_name}' ya existe. Eliminando...")
            shutil.rmtree(new_folder_path)
        
        shutil.copytree(source_folder, new_folder_path)
        
        # Modificar el archivo input/pump
        pump_file = new_folder_path / "input" / "pump"
        
        if pump_file.exists():
            with open(pump_file, 'w') as f:
                f.write(f"{pump_value}\n")
            print(f"  ✓ Archivo 'input/pump' actualizado con valor: {pump_value}")
        else:
            print(f"  ⚠️  Advertencia: No se encontró el archivo 'input/pump' en {new_folder_name}")
        
        print()
    
    print("¡Proceso completado!")
    print(f"\nCarpetas creadas:")
    for pump_value in pump_values:
        print(f"  - pump={pump_value}")


if __name__ == "__main__":
    main()
