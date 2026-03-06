import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import TS_calculator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import matplotlib.pyplot as plt
import os

# --- CONFIGURATION INITIALE ---
date = 'fevrier'


freq = "70kHz"


if freq == "70kHz":
    channel = 1
    if date == 'octobre':
        calib_file = "EK80_Panopee/Calib_octobre/CalibrationDataFile-D20251001-T143810_70kHz_OK.xml"
        output_base_dir = "Analyse_Massive_Cibles_70kHz_octobre"
        input_root = "EK80_Panopee/Octobre_raws/0910"

    else : 
        calib_file = "EK80_Panopee/Calib_fevrier/0402_70kHz/CalibrationDataFile-D20260204-T134413-70kHz.xml"
        output_base_dir = "Analyse_Massive_Cibles_70kHz_fevrier"
        input_root = "EK80_Panopee/Fevrier_raws/0502"
else :
    channel = 0
    if date == 'octobre':
        calib_file = "EK80_Panopee/Calib_octobre/CalibrationDataFile-D20251001-T140822_200kHz_OK.xml"
        output_base_dir = "Analyse_Massive_Cibles_200kHz_octobre"
        input_root = "EK80_Panopee/Octobre_raws"
    else : 
        calib_file = "EK80_Panopee/Calib_fevrier/0302_200kHz/CalibrationDataFile-D20260203-T101032-200kHz.xml"
        output_base_dir = "Analyse_Massive_Cibles_200kHz_fevrier"
        input_root = "EK80_Panopee/Fevrier_raws"

# 1. Initialisation du dictionnaire GLOBAL (hors de la fonction)
trajectoires = {}
file_ids = []


def process_file(file_path): 
    rel_path = os.path.relpath(file_path, input_root)
    file_id = rel_path.replace(os.sep, "_").replace(".raw", "")
    
    try:
        proc = TS_calculator.EK80Processor(calib_file=calib_file, channel=channel)
        proc.load_raw(file_path)

        x_raw = proc.x_boat_gps
        y_raw = proc.y_boat_gps

        # FILTRAGE : On ne garde que les points dans un rectangle cohérent autour du lac
        # Cela élimine les sauts vers (0,0) qui créent les lignes droites
        mask = (x_raw > 200000) & (y_raw > 6000000) 
        
        x_clean = x_raw[mask]
        y_clean = y_raw[mask]

        if len(x_clean) > 0:
            trajectoires[file_id] = {
                'x': x_clean,
                'y': y_clean
            }
            file_ids.append(file_id)
            print(f"Position extraite et filtrée pour : {file_id}")
        else:
            print(f"⚠️ Aucun point valide dans la zone pour : {file_id}")

    except Exception as e:
        print(f"Erreur lors du traitement de {rel_path} : {e}")


# --- EXÉCUTION ---

# Recherche de tous les fichiers .raw
raw_files = sorted(glob.glob(os.path.join(input_root, "**", "*.raw"), recursive=True))

# Boucle de traitement (remplit le dictionnaire)
for f in raw_files:
    process_file(f)

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import matplotlib.pyplot as plt
import os

def tracer_levee_cartopy(trajectoires, output_base_dir):
    # 1. Définition de la projection Lambert-93 pour Cartopy
    # Le Lambert-93 est une projection "Lambert Conformal" spécifique à la France
    projection_l93 = ccrs.LambertConformal(
        central_longitude=3.0, 
        central_latitude=46.5,
        standard_parallels=(44.0, 49.0),
        false_easting=700000,
        false_northing=6600000
    )

    fig = plt.figure(figsize=(12, 10))
    # On crée un axe avec la projection Lambert-93
    ax = plt.axes(projection=projection_l93)

    # 2. Ajout d'un fond de carte (OpenStreetMap)
    # On choisit un niveau de zoom élevé (14 ou 15) pour bien voir le lac
    osm_tiles = cimgt.OSM()
    ax.add_image(osm_tiles, 14, alpha=0.6)

    # 3. Tracé des trajectoires
    # Comme tes données sont déjà en Lambert-93 (mètres), on précise transform=projection_l93
    colors = plt.cm.viridis(np.linspace(0, 1, len(trajectoires)))
    
    for (file_id, coords), color in zip(trajectoires.items(), colors):
        ax.plot(coords['x'], coords['y'], color=color, lw=1.5, 
                transform=projection_l93, label=file_id.split('-')[-1])
        print(np.mean(coords['x']))
        

    # 4. Définition de la zone de vue (Zoom sur Guerlédan)
    # Coordonnées Lambert-93 min/max
    plt.xlim(252500,254000)

    x_mid = [253327.79509856337, 252940.27602594183, 252248.87139314972]
    y_mid = [6805617.805236985, 6805526.134501897,  6805765.435748426 ]
    zones = ['Zone 1', 'Zone 2', 'Zone 3']

    for i in range(len(x_mid)):
        zone_label = zones[i]
        ax.text(x_mid[i] + 50, y_mid[i] + 50, zone_label, 
        transform=projection_l93, fontsize=8, fontweight='bold',
        bbox=dict(facecolor='white', alpha=0.6, edgecolor='black', boxstyle='round'))


    ax.set_extent([251500, 254000, 6804800, 6806200], crs=projection_l93)

    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha=0.3)
    plt.title(f"Plan de levé - Lac de Guerlédan", fontweight='bold')
    
    # Légende (si pas trop de fichiers)
    if len(trajectoires) < 10:
        plt.legend(loc='lower right', fontsize='small')

    save_path = os.path.join(output_base_dir, "Carte_Cartopy_Guerledan.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

# Appel de la fonction
tracer_levee_cartopy(trajectoires, output_base_dir)



