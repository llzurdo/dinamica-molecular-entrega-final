import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del archivo
data = np.loadtxt('output/perfil_velocidad.dat', comments='#')

# Separar las columnas
z_centro = data[:, 0]
vx_promedio = data[:, 1]

# Crear el gráfico
plt.figure(figsize=(10, 6))
plt.plot(z_centro, vx_promedio, marker='o', linestyle='-', color='b')
plt.title('Perfil de Velocidad en Función de Z')
plt.xlabel('Z Centro')
plt.ylabel('Vx Promedio')
plt.grid(True)
plt.savefig("output/perfil_velocidad.pdf", dpi=300)
