import matplotlib.pyplot as plt
import TS_calculator
import os
import numpy as np

calib_file = "EK80_Panopee/Calib_fevrier/0402_70kHz/CalibrationDataFile-D20260204-T134413-70kHz.xml"
channel = 1  # 0 pour 200kHz, 1 pour 70kHz
date = "fevrier"
output_dir = "EK80_Panopee/calib_plots"

def plot_gain_calib(file_path, channel, date, output_dir):
    freq_label = "70kHz" if channel == 1 else "200kHz"
    

    proc = TS_calculator.EK80Processor(calib_file=file_path, channel=channel)

    freqs = np.array(proc.frequencies)

    plt.figure(figsize=(10, 6))
    plt.plot(freqs, proc.gain, color='firebrick', lw=2, label='Gain de Calibration')

    plt.title(f"Gain du Transducteur - {freq_label} ({date})", fontsize=14)
    plt.xlabel("Fréquence (kHz)", fontsize=12)
    plt.ylabel("Gain (dB)", fontsize=12)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend()

    # Création propre du dossier
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Dossier créé : {os.path.abspath(output_dir)}")

    filename = f"Gain_Calibration_{date}_{freq_label}.png"
    save_path = os.path.join(output_dir, filename)
    
    plt.savefig(save_path, dpi=150)
    plt.close()

plot_gain_calib(calib_file,channel,date,output_dir)
