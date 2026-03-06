import numpy as np
import matplotlib.pyplot as plt

def simulate_krm_sigma_unique(L, tilt_deg, freq_range):
    c = 1490.0
    L_v, a_v = L * 0.30, L * 0.025 
    theta = np.radians(tilt_deg)
    k = 2 * np.pi * freq_range / c

    k_ref = 2 * np.pi * 70000 / c
    
    X = k * L_v * np.sin(theta)
    directivity = np.sinc(X / np.pi)
    R_coef = -0.99
    
    sigma = ((k_ref * a_v * L_v**2) / (4 * np.pi)) * (directivity**2) * (np.cos(theta)**2) * np.abs(R_coef)**2
    return 10 * np.log10(sigma + 1e-12)

def simulate_krm_sigma_complexe(L, tilt_deg, freq_range):
    c = 1490.0
    R = -0.99 
    theta = np.radians(tilt_deg)
    k = 2 * np.pi * freq_range / c
    k_ref = 2 * np.pi * 70000 / c
    
    dist = L * 0.5 
    L1, a1, z1 = L * 0.20, L * 0.030, 0.010
    L2, a2, z2 = L * 0.15, L * 0.020, 0.0
    
    def get_p(Li, ai, x, z):
        X = k * Li * np.sin(theta)
        # On utilise sqrt(k_ref) pour que l'amplitude soit stable en fréquence
        amp = np.sqrt((k_ref * ai * Li**2) / (4 * np.pi)) * R * np.sinc(X/np.pi) * np.cos(theta)
        # La PHASE reste dépendante de k (car c'est elle qui crée les interférences)
        phase = np.exp(2j * k * (x * np.sin(theta) + z * np.cos(theta)))
        return amp * phase

    p_total = get_p(L1, a1, -dist/2, z1) + get_p(L2, a2, dist/2, z2)
    return 10 * np.log10(np.abs(p_total)**2 + 1e-12)

# --- Simulation ---
L_fish = 0.2

freqs = np.linspace(45000, 95000, 1000)
angles_etude = np.arange(0, -5, -1)  # Pas de 2 degrés

# --- Affichage ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 10), sharex=True)
colors = plt.cm.plasma(np.linspace(0, 1, len(angles_etude)))

# Subplot 1 : Vessie Unique
for i, angle in enumerate(angles_etude):
    ts_spectrum = simulate_krm_sigma_unique(L_fish, angle, freqs)
    ax1.plot(freqs/1000, ts_spectrum, label=f'Angle: {angle+90}°', color=colors[i], lw=1.5)
ax1.set_ylabel('Target Strength (dB)',fontsize=24)
#ax1.set_title('Modèle KRM : Vessie à chambre unique')
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=24)

# Subplot 2 : Double Chambre
for i, angle in enumerate(angles_etude):
    ts_spectrum = simulate_krm_sigma_complexe(L_fish, angle, freqs)
    ax2.plot(freqs/1000, ts_spectrum, label=f'Angle: {angle+90}°', color=colors[i], lw=1.5)
ax2.set_xlabel('Fréquence (kHz)',fontsize=24)
ax2.set_ylabel('Target Strength (dB)',fontsize=24)
#ax2.set_title('Modèle KRM : Vessie à double chambre (Physostome)')
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=24)

plt.tight_layout()
plt.savefig('KRM')
plt.close()

import numpy as np
import matplotlib.pyplot as plt

# --- Fonctions de base modifiées pour retourner SIGMA (linéaire) ---

def simulate_krm_sigma_unique(L, tilt_deg, freq_range):
    c = 1490.0
    L_v, a_v = L * 0.30, L * 0.025 
    theta = np.radians(tilt_deg)
    k = 2 * np.pi * freq_range / c

    k_ref = 2 * np.pi * 70000 / c
    
    X = k * L_v * np.sin(theta)
    directivity = np.sinc(X / np.pi)
    R_coef = -0.99
    
    sigma = ((k_ref * a_v * L_v**2) / (4 * np.pi)) * (directivity**2) * (np.cos(theta)**2) * np.abs(R_coef)**2
    return sigma

def simulate_krm_sigma_complexe(L, tilt_deg, freq_range):
    c = 1490.0
    R = -0.99 
    theta = np.radians(tilt_deg)
    k = 2 * np.pi * freq_range / c
    k_ref = 2 * np.pi * 70000 / c
    
    dist = L * 0.5 
    L1, a1, z1 = L * 0.20, L * 0.030, 0.010
    L2, a2, z2 = L * 0.15, L * 0.020, 0.0
    
    def get_p(Li, ai, x, z):
        X = k * Li * np.sin(theta)
        # On utilise sqrt(k_ref) pour que l'amplitude soit stable en fréquence
        amp = np.sqrt((k_ref * ai * Li**2) / (4 * np.pi)) * R * np.sinc(X/np.pi) * np.cos(theta)
        # La PHASE reste dépendante de k (car c'est elle qui crée les interférences)
        phase = np.exp(2j * k * (x * np.sin(theta) + z * np.cos(theta)))
        return amp * phase

    p_total = get_p(L1, a1, -dist/2, z1) + get_p(L2, a2, dist/2, z2)
    return np.abs(p_total)**2

def get_averaged_heatmap(model_func, L_start, L_end, steps, angles, freqs):
    L_range = np.linspace(L_start, L_end, steps)
    sum_sigma = np.zeros((len(angles), len(freqs)))
    for L in L_range:
        for i, ang in enumerate(angles):
            # Calcul en linéaire
            sum_sigma[i, :] += model_func(L, ang, freqs)
    # Moyenne linéaire puis passage en dB
    avg_sigma = sum_sigma / steps
    return 10 * np.log10(avg_sigma + 1e-12)

# --- Paramètres de simulation ---
freqs = np.linspace(45000, 95000, 400)
angles = np.linspace(-10, 10, 150) # +/- 10° autour de 90°

# Calcul des heatmaps moyennées (3cm à 20cm)
data_unique = get_averaged_heatmap(simulate_krm_sigma_unique, 0.10, 0.15, 30, angles, freqs)
data_complexe = get_averaged_heatmap(simulate_krm_sigma_complexe, 0.10, 0.15, 30, angles, freqs)

# --- Affichage ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

for ax, data, title in zip([ax1, ax2], [data_unique, data_complexe], 
                           ['Vessie Unique (Moy. 3-20cm)', 'Double Chambre (Moy. 3-20cm)']):
    im = ax.imshow(data, aspect='auto', origin='lower',
                   extent=[freqs[0]/1000, freqs[-1]/1000, 80, 100], cmap='viridis')
    #ax.set_title(title, fontsize=16)
    ax.set_xlabel('Fréquence (kHz)',fontsize=24)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.legend(fontsize=24)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('TS Moyen (dB re 1m²)', fontsize=24, labelpad=15)
cbar.ax.tick_params(labelsize=20)

ax1.set_ylabel('Angle d\'inclinaison (°)',fontsize=24)

plt.tight_layout()
plt.savefig('KRM_heatmap_moyenne.png')
