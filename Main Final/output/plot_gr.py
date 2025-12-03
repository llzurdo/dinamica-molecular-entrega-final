import numpy as np
import matplotlib.pyplot as plt

# Cargar gr.dat ignorando la línea de comentario
data = np.loadtxt("gr.dat", comments="#")

# La estructura del archivo es:
#   r   g(r)   counts
r = data[:, 0]
g = data[:, 1]

plt.figure(figsize=(8, 5))
plt.plot(r, g, linewidth=2)

plt.xlabel("r", fontsize=14)
plt.ylabel("g(r)", fontsize=14)
plt.title("Función de distribución radial g(r)", fontsize=16)

plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

# Guardar la figura
plt.savefig("gr.png", dpi=150)

# Mostrar en pantalla
plt.show()
