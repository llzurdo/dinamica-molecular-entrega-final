import subprocess
import os
import sys

def convert_to_mp4(input_file, output_file=None, quality="medium", overwrite=False):
    """
    Convierte un video a formato MP4 usando FFmpeg.
    
    Args:
        input_file (str): Ruta del archivo de video de entrada
        output_file (str): Ruta del archivo MP4 de salida (opcional)
        quality (str): Calidad del video ("low", "medium", "high")
        overwrite (bool): Si True, sobrescribe el archivo de salida si existe
    
    Returns:
        bool: True si la conversión fue exitosa, False en caso contrario
    """
    
    # Verificar que el archivo de entrada existe
    if not os.path.exists(input_file):
        print(f"❌ Error: El archivo '{input_file}' no existe.")
        return False
    
    # Generar nombre de salida si no se proporciona
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_converted.mp4"
    
    # Verificar si el archivo de salida ya existe
    if os.path.exists(output_file) and not overwrite:
        print(f"❌ Error: El archivo '{output_file}' ya existe. Use overwrite=True para sobrescribir.")
        return False
    
    # Configuraciones de calidad
    quality_settings = {
        "low": ["-crf", "28", "-preset", "fast"],
        "medium": ["-crf", "23", "-preset", "medium"],
        "high": ["-crf", "18", "-preset", "slow"]
    }
    
    if quality not in quality_settings:
        print(f"❌ Error: Calidad '{quality}' no válida. Use: low, medium, high")
        return False
    
    # Construir comando FFmpeg
    cmd = [
        "ffmpeg",
        "-i", input_file,
        "-c:v", "libx264",           # Codec de video H.264
        "-c:a", "aac",               # Codec de audio AAC
        "-movflags", "+faststart",   # Optimizar para streaming web
    ] + quality_settings[quality] + [
        "-y" if overwrite else "-n", # Sobrescribir o no
        output_file
    ]
    
    try:
        print(f"🔄 Convirtiendo '{input_file}' a MP4...")
        print(f"📁 Archivo de salida: '{output_file}'")
        print(f"⚙️ Calidad: {quality}")
        
        # Ejecutar FFmpeg
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"✅ Conversión exitosa!")
            print(f"📂 Archivo guardado en: {os.path.abspath(output_file)}")
            return True
        else:
            print(f"❌ Error en la conversión:")
            print(f"Error: {result.stderr}")
            return False
            
    except FileNotFoundError:
        print("❌ Error: FFmpeg no está instalado o no está en el PATH.")
        print("Para instalar FFmpeg en Linux: sudo apt install ffmpeg")
        return False
    except Exception as e:
        print(f"❌ Error inesperado: {e}")
        return False

def batch_convert_to_mp4(input_dir, output_dir=None, quality="medium", overwrite=False):
    """
    Convierte múltiples videos a MP4 en lote.
    
    Args:
        input_dir (str): Directorio con videos de entrada
        output_dir (str): Directorio de salida (opcional, usa el mismo si no se especifica)
        quality (str): Calidad del video
        overwrite (bool): Sobrescribir archivos existentes
    """
    
    if not os.path.exists(input_dir):
        print(f"❌ Error: El directorio '{input_dir}' no existe.")
        return
    
    if output_dir is None:
        output_dir = input_dir
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Extensiones de video comunes
    video_extensions = {'.avi', '.mov', '.mkv', '.flv', '.wmv', '.webm', '.m4v', '.3gp'}
    
    video_files = [f for f in os.listdir(input_dir) 
                   if os.path.splitext(f)[1].lower() in video_extensions]
    
    if not video_files:
        print(f"❌ No se encontraron archivos de video en '{input_dir}'")
        return
    
    print(f"🎬 Encontrados {len(video_files)} archivos de video")
    
    successful = 0
    for video_file in video_files:
        input_path = os.path.join(input_dir, video_file)
        output_filename = os.path.splitext(video_file)[0] + ".mp4"
        output_path = os.path.join(output_dir, output_filename)
        
        print(f"\n--- Procesando {video_file} ---")
        if convert_to_mp4(input_path, output_path, quality, overwrite):
            successful += 1
    
    print(f"\n🎯 Conversión completada: {successful}/{len(video_files)} archivos exitosos")

# ===== EJEMPLO DE USO =====
if __name__ == "__main__":
    # Ejemplo 1: Convertir un solo archivo
    # convert_to_mp4("video_original.avi", "video_convertido.mp4", quality="high")
    
    # Ejemplo 2: Conversión en lote
    # batch_convert_to_mp4("videos_input/", "videos_output/", quality="medium")
    
    # Ejemplo 3: Convertir el video de tu simulación
    if os.path.exists("output/animation.mp4"):
        print("El archivo ya está en formato MP4")
    else:
        print("No se encontró el archivo de animación para convertir")