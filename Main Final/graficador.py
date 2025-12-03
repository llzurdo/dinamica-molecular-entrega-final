from ovito.io import import_file
from ovito.vis import Viewport, CoordinateTripodOverlay, ViewportOverlay, TextLabelOverlay
import OpenGL.GL as gl
import os
from utils.convertidor_a_mp4 import convert_to_mp4
import sys
import numpy as np
import matplotlib.pyplot as plt

'''
Parte que grafica la energia
'''

fn = "output/energia.dat"

with open("input/deltaT", "r") as f:
    deltaT = float(f.read().strip())

with open("input/deltaTMinimizacion", "r") as f:
    deltaTMinimizacion = float(f.read().strip())
    
with open("input/muestreo", "r") as f:
    muestreo = int(f.read().strip())

with open("input/muestreoMinimizacion", "r") as f:
    muestreoMinimizacion = int(f.read().strip())

with open("input/tMinimizacion", "r") as f:
    tMinimizacion = int(int(f.read().strip()) / muestreoMinimizacion)

try:
    data = np.loadtxt(fn)
except Exception as e:
    print(f"Error al leer '{fn}': {e}", file=sys.stderr)
    sys.exit(1)

# detectar columnas: si hay >=2 columnas, usar primera como x y segunda como y; si sólo hay 1, usar índice como x
if data.ndim == 1:
    y_minimizacion = data[:tMinimizacion]
    x_minimizacion = np.arange(len(y_minimizacion)) * deltaTMinimizacion * muestreoMinimizacion
    y              = data[tMinimizacion:]
    x              = np.arange(len(y)) * deltaT * muestreo + max(x_minimizacion)

# Para concatenar vectores 
x_final = np.concatenate([x_minimizacion, x])
y_final = np.concatenate([y_minimizacion, y])

# Para gráficos
mean    = np.mean(y)
std     = np.std(y)
ymin    = -np.abs(mean*1.5)-np.abs(std);
ymax    = np.abs(mean*1.5)+np.abs(std);

plt.figure(figsize=(8, 4.5), dpi = 300)
plt.plot(x_final, y_final, '-', markersize=3)
plt.axhline(y=mean, color='red', linestyle='--', linewidth=1, label=f'Mean = {mean:.2f}, Std: {std:.3f}')
plt.xlabel('Tiempo' if data.ndim == 1 or data.shape[1] == 1 else 'columna 1')
plt.ylabel('Energia' if data.ndim > 1 and data.shape[1] > 1 else 'Energía')
plt.title('Energía vs Tiempo' + fn)
plt.ylim(ymin, ymax)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("output/energia.pdf")


'''
Parte que grafica la presion
'''

fn = "output/presion.dat"

try:
    data = np.loadtxt(fn)
except Exception as e:
    print(f"Error al leer '{fn}': {e}", file=sys.stderr)
    sys.exit(1)

# detectar columnas: si hay >=2 columnas, usar primera como x y segunda como y; si sólo hay 1, usar índice como x
if data.ndim == 1:
    y = data;
    x = np.arange(len(y)) * deltaT * muestreo

# Para gráficos
mean    = np.mean(y)
std     = np.std(y)
ymin    = -np.abs(mean*1.5)-np.abs(std);
ymax    = np.abs(mean*1.5)+np.abs(std);

plt.figure(figsize=(8, 4.5), dpi = 300)
plt.plot(x, y, '-', markersize=3)
plt.axhline(y=mean, color='red', linestyle='--', linewidth=1, label=f'Mean = {mean:.5f}, Std: {std:.5f}')
plt.xlabel('Tiempo' if data.ndim == 1 or data.shape[1] == 1 else 'columna 1')
plt.ylabel('Presion' if data.ndim > 1 and data.shape[1] > 1 else 'Energía')
plt.title('Presion vs Tiempo' + fn)
plt.ylim(ymin, ymax)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("output/presion.pdf")


'''
Parte que grafica la Temperatura
'''

fn = "output/temperatura.dat"

try:
    data = np.loadtxt(fn)
except Exception as e:
    print(f"Error al leer '{fn}': {e}", file=sys.stderr)
    sys.exit(1)

# detectar columnas: si hay >=2 columnas, usar primera como x y segunda como y; si sólo hay 1, usar índice como x
if data.ndim == 1:
    y = data;
    x = np.arange(len(y)) * deltaT * muestreo

# Para gráficos
mean    = np.mean(y)
std     = np.std(y)
ymin    = -np.abs(mean*1.5)-np.abs(std);
ymax    = np.abs(mean*1.5)+np.abs(std);

plt.figure(figsize=(8, 4.5), dpi = 300)
plt.plot(x, y, '-', markersize=3)
plt.axhline(y=mean, color='red', linestyle='--', linewidth=1, label=f'Mean = {mean:.4f}, Std: {std:.4f}')
plt.xlabel('Tiempo' if data.ndim == 1 or data.shape[1] == 1 else 'columna 1')
plt.ylabel('Temperatura' if data.ndim > 1 and data.shape[1] > 1 else 'Energía')
plt.title('Temperatura vs Tiempo' + fn)
plt.ylim(ymin, ymax)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("output/temperatura.pdf")

'''

# Graficamos el video
# ===== CONFIGURACIÓN =====
input_file      = "output/posicionT.xyz"
output_dir      = "output"
image_size      = (800, 600)
with open("input/sigma", "r") as f:
    particle_radius = float(f.read().strip())/2
    print(f"Radio de partícula leído desde archivo: {particle_radius}")
scale_length    = 10.0  # 🔹 Longitud de la barra de escala en unidades de simulación
bar_color       = (0.0, 0.0, 0.0)  # negro

# ===== CARGA DEL ARCHIVO =====
print(f"Cargando archivo: {input_file}")
pipeline = import_file(input_file)

num_frames = pipeline.source.num_frames
print(f"Frames detectados: {num_frames}")

# ===== CONFIGURAR VISUALIZACIÓN =====
pipeline.add_to_scene()

particles_vis = pipeline.source.data.particles.vis
particles_vis.radius = particle_radius

# ===== CREAR CARPETA =====
os.makedirs(output_dir, exist_ok=True)

# ===== CONFIGURAR VIEWPORT =====
vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.zoom_all()

# Alejar la cámara
pos           = vp.camera_pos
vp.camera_pos = (pos[0] * 0.8 + 0.1, pos[1] * 0.8 + 0.1, pos[2] * 0.8 + 0.1)

# ===== AGREGAR TRÍPODE =====
tripod            = CoordinateTripodOverlay()
tripod.size       = 0.03
tripod.line_width = 0.5
vp.overlays.append(tripod)


output_path = os.path.join(output_dir, f"animation.mp4")

#print(f"Renderizando video → {output_path}")
#vp.render_anim(size=(800,600), filename=output_path, fps=20)


#print("\n✅ Renderización completa.")
#print(f"Video guardado en: {os.path.abspath(output_dir)}")
'''
# Convertimos el video
'''
convert_to_mp4(input_file="output/animation.mp4",
               overwrite=True,
               output_file="output/animation2.mp4",
               quality="medium")
'''
