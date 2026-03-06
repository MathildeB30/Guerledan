import matplotlib.pyplot as plt
import gc
import TS_calculator
import numpy as np
import pandas as pd
import os
import glob
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import interp1d
import seaborn as sns
from scipy.interpolate import RegularGridInterpolator

date = 'octobre'



freq="70kHz"



if freq == "70kHz":
    channel = 1
    if date == 'octobre':
        calib_file = "EK80_Panopee/Calib_octobre/CalibrationDataFile-D20251001-T143810_70kHz_OK.xml"
        output_base_dir = "Analyse_Massive_Cibles_70kHz_octobre"
        input_root = "EK80_Panopee/Octobre_raws"
    else : 
        calib_file = "EK80_Panopee/Calib_fevrier/0402_70kHz/CalibrationDataFile-D20260204-T134413-70kHz.xml"
        output_base_dir = "Analyse_Massive_Cibles_70kHz_fevrier"
        input_root = "EK80_Panopee/Fevrier_raws"
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


colors_list = plt.cm.tab10.colors 

TS_MIN = -80
PROM_MIN = 10
WIDTH_RANGE = (8, 100)
DIST_MAX_M = 2.5
TIME_BETWEEN_PINGS = 0.5

raw_files = sorted(glob.glob(os.path.join(input_root, "**", "*.raw"), recursive=True))
print(f"Trouvé {len(raw_files)} fichiers .raw.")

all_targets_features = []
resultats_globaux = []


def calc_mean_ts(group):
    ts_spectra = [d['ts'] for d in group]
    ts_linear = 10**(np.array(ts_spectra) / 10.0)
    mean_linear = np.mean(ts_linear, axis=0)
    return 10 * np.log10(mean_linear)

def plot_echogram(sp_matrix, r_n, output_path,title):
    fig, ax1 = plt.subplots(figsize=(10, 8))
    extent = [0, sp_matrix.shape[0], r_n[-1], r_n[0]]
    im1 = ax1.imshow(sp_matrix.T, aspect='auto', cmap='jet', vmin=-100, vmax=-20, extent=extent)
    ax1.set_xlabel('Numéro de Ping', fontsize=24, labelpad=10)
    ax1.set_ylabel('Profondeur (m)', fontsize=24, labelpad=10)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    cbar = fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label('Target Strength (dB re 1m²)', fontsize=24, labelpad=15)
    cbar.ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

def plot_spectra(target,f_m, output_path,incident_angles):
    fig, ax = plt.subplots(figsize=(10, 8))

    for i, t in enumerate(target):
        color = colors_list[i % 10]
        ax.plot(f_m/1000, t['ts'], color=color, alpha=0.5, label=f"P{t['ping']} ({incident_angles[i]:.1f}°)")
    ax.set_title("Spectres TSf de la cible",fontsize=24)
    ax.set_xlabel("Fréquence (kHz)",fontsize=24)
    ax.set_ylabel("TS (dB)",fontsize=24)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.legend(fontsize=18, ncol=2)
    ax.grid(True, alpha=0.2)
    plt.savefig(output_path, dpi=150)
    plt.close()


def process_file(file_path): 
    rel_path = os.path.relpath(file_path, input_root)
    file_id = rel_path.replace(os.sep, "_").replace(".raw", "")
    
    print(f"\n>>> Traitement de : {rel_path}")
    file_output_dir = os.path.join(output_base_dir, file_id)
    if not os.path.exists(file_output_dir): os.makedirs(file_output_dir)


    try:
        proc = TS_calculator.EK80Processor(calib_file=calib_file, channel=channel)
        proc.load_raw(file_path)

        f_s_dec = proc.calcDecmiatedSamplingRate()
        print(f_s_dec,proc.f_s)
        tau = proc.tau
    

        
        samples_pulse = tau * f_s_dec
        WINDOW_FFT = int(2**np.ceil(np.log2(samples_pulse * 2)))
            
        all_pings_data = [] 
        all_peaks_data = []
        SP_all = []
        SP_all_full = []
        indices_fond = []
        
        for ping in range(proc.nb_pings):
            y_rx_nu = (proc.y_rx_nu[channel,ping,:,:]).T
            y_rx_nu = proc._trim_nans(y_rx_nu, axis=1)
            sampleCount = y_rx_nu.shape[1]
            r_n = proc.calcRange(sampleCount)
            y_pc_nu = proc.calcPulseCompressedSignals(y_rx_nu)

            y_pc_halves = proc.calcTransducerHalves(y_pc_nu)
            y_pc = proc.calcAverageSignal(y_pc_nu)
            p_rx_n = proc.calcPower(y_pc)



            theta, phi = proc.calcAngles(y_pc_halves)

            limit = 4253 if channel == 1 else 8507
            sp = proc.calcSp(p_rx_n, theta, phi, r_n)[:limit]
            r_n = r_n[:limit]

            start_search = 150 
            idx_hard_bottom = np.argmax(sp[start_search:]) + start_search
            threshold_bottom = -100
            search_zone_up = sp[start_search:idx_hard_bottom]
            idx_below_thresh = np.where(search_zone_up < threshold_bottom)[0]
            
            if len(idx_below_thresh) > 0:
                idx_fond = idx_below_thresh[-1] + start_search
            else:
                idx_fond = idx_hard_bottom - 5

            indices_fond.append(idx_fond)
            sp_clean = sp.copy()
            sp_clean[idx_fond:] = -150 
            sp_clean[:100] = -150
            SP_all_full.append(sp)
            SP_all.append(sp_clean)

            peaks, _ = find_peaks(sp_clean, height=TS_MIN, prominence=PROM_MIN, width=WIDTH_RANGE)
            
            if len(peaks) > 0:
                widths, _, left_ips, right_ips = peak_widths(sp_clean, peaks, rel_height=0.5)
                for i in range(len(peaks)):
                    all_peaks_data.append({
                        'ping_index': ping,
                        'max_idx': peaks[i],
                        'sp' : sp_clean[peaks[i]],      
                        'start_idx': int(np.round(left_ips[i])),
                        'end_idx': int(np.round(right_ips[i]))
                    })

            all_pings_data.append({'y_pc_n': y_pc, 'theta_n': theta, 'phi_n': phi, 'r_n': r_n})

        SP_all = np.array(SP_all)
        SP_all_full = np.array(SP_all_full)

        plot_echogram(SP_all_full, r_n, os.path.join(file_output_dir, "echogram_full.png"), "Echogramme Complet")
        plot_echogram(SP_all, r_n, os.path.join(file_output_dir, "echogram_clean.png"), "Echogramme Nettoyé (Fond Masqué)")
            
        df_peaks = pd.DataFrame(all_peaks_data)
        groups = []

        for i in range(len(all_peaks_data)):
            ping = int(df_peaks['ping_index'][i])
            data = all_pings_data[ping]
            max_idx = df_peaks['max_idx'][i]
            r_t = r_n[max_idx]

            start_win = int(max(0, max_idx - WINDOW_FFT // 2))
            end_win = int(start_win + WINDOW_FFT)

            if end_win > len(data['y_pc_n']):
                end_win = len(data['y_pc_n'])
                start_win = max(0, end_win - WINDOW_FFT)
                

            y_pc_fixed = data['y_pc_n'][start_win:end_win]
            
            theta_t, phi_t = data['theta_n'][max_idx], data['phi_n'][max_idx]
            x_loc = r_t * np.sin(np.deg2rad(theta_t))
            y_loc = r_t * np.sin(np.deg2rad(phi_t))
            z_real = r_t * np.cos(np.sqrt(np.deg2rad(theta_t)**2 + np.deg2rad(phi_t)**2))
            x_abs = proc.x_boat_gps[ping] + x_loc
            y_abs = proc.y_boat_gps[ping] + y_loc
            pos_rel = np.array([x_loc, y_loc, z_real])

            y_mf_auto_red_n = proc.alignAuto(y_pc_fixed)
            _, _, Y_tilde = proc.calcDFTforTS(y_pc_fixed, y_mf_auto_red_n, f_s_dec)
            TS_m = proc.calcTSf(proc.calcPowerFreqTS(Y_tilde), r_t, theta_t, phi_t)

            current_det = {
                'ping': ping, 'sp' : df_peaks['sp'][i], 'ts': TS_m, 'x_abs': x_abs, 'y_abs': y_abs, 
                'z': z_real, 'pos_rel': pos_rel, 'theta_n': theta_t, 'phi_n': phi_t
            }

            assigned = False
            for group in groups :
                last_det = group[-1]
                dz = np.abs(current_det['z'] - last_det['z'])
                dist = np.sqrt((current_det['x_abs'] - last_det['x_abs'])**2 + (current_det['y_abs'] - last_det['y_abs'])**2 + (current_det['z'] - last_det['z'])**2)
                time_gap = current_det['ping'] - last_det['ping']

                if dist <= DIST_MAX_M and time_gap <= 2 and dz < 0.1:
                    group.append(current_det)
                    assigned = True
                    break
                
            if not assigned:
                groups.append([current_det])
        print(f"Sauvegardé : {len(groups)} cibles")
        for idx, target in enumerate(groups):
            if len(target) >= 3: 

                incident_angles = []
                incident_angles.append(np.nan)
                total_dist = 0
                for i in range (len(target)-1):
                    v_fish = np.array([
                        target[i+1]['x_abs'] - target[i]['x_abs'],
                        target[i+1]['y_abs'] - target[i]['y_abs'],
                        target[i+1]['z'] - target[i]['z']
                    ])
                    total_dist += np.linalg.norm(v_fish)
            
                    v_beam = target[i+1]['pos_rel'] 
                    cross_prod = np.cross(v_beam, v_fish)
                    norm_cross = np.linalg.norm(cross_prod)
                    norm_X = np.linalg.norm(v_beam)
                    norm_Z = np.linalg.norm(v_fish)
                    sin_theta = norm_cross / (norm_X * norm_Z)
                    theta_rad = np.arcsin(np.clip(sin_theta, -1.0, 1.0))
                    angle_deg = np.rad2deg(theta_rad)
                    incident_angles.append(angle_deg)
                
                total_time = (target[-1]['ping'] - target[0]['ping']) * TIME_BETWEEN_PINGS
                velocity_mps = total_dist / total_time if total_time > 0 else 0

                mean_ts = calc_mean_ts(target)
                
                if np.nanmax(mean_ts) < -60:
                    continue

                all_ts_specs = np.array([t['ts'] for t in target]) 
                power_per_ping = np.nanmean(10**(all_ts_specs / 10.0), axis=1)
                best_ping_idx = np.argmax(power_per_ping)
                ts_at_peak_db = 10 * np.log10(power_per_ping[best_ping_idx])
                tilt_angle_peak = incident_angles[best_ping_idx]
                
                
                ts_diff = [np.diff(target['ts']) for target in target]
                spectral_rugosity = [np.std(ts_d) for ts_d in ts_diff]
                nb_high_rugosity = sum(1 for x in spectral_rugosity if x > 3)
                majorite_seuil = len(spectral_rugosity) / 2

                if nb_high_rugosity > majorite_seuil:
                    save_path = os.path.join(file_output_dir, f"Target_{idx+1:03d}_rejected_fond.png")
                    is_rejected = True
                else : 
                    # Conservé comme poisson potentiel
                    save_path = os.path.join(file_output_dir, f"Target_{idx+1:03d}_fish.png")
                    is_rejected = False


                # Préparation des données pour les graphiques
                xs = np.array([d['x_abs'] for d in target])
                ys = np.array([d['y_abs'] for d in target])
                zs = np.array([d['z'] for d in target])
                pings = np.array([d['ping'] for d in target])
                xb_target = proc.x_boat_gps[pings]
                yb_target = proc.y_boat_gps[pings]
                dist_fish_2d = np.sqrt((xs+xb_target)**2 + (ys+yb_target)**2)

                n_angle_std = 30
                n_freq_std = len(proc.f_m)
                std_angles = np.linspace(70, 100, n_angle_std)
                std_freqs = proc.f_m/1000

                raw_angles = np.array(incident_angles)
                raw_specs = np.array([t['ts'] for t in target])
                mask = ~np.isnan(raw_angles)
                x_clean = raw_angles[mask]
                y_clean = raw_specs[mask]
                sort_idx = np.argsort(x_clean)
                
                x_sorted = x_clean[sort_idx]
                y_sorted = y_clean[sort_idx]

                f_spat = interp1d(x_sorted, y_sorted, axis=0, kind='linear', 
                                bounds_error=False, 
                                fill_value="extrapolate")

                
                matrix_intermediate = f_spat(std_angles)

                f_freq = interp1d(proc.f_m / 1000.0, matrix_intermediate, axis=1, 
                                kind='linear', 
                                bounds_error=False, 
                                fill_value="extrapolate") # L'extrapol. freq est moins risquée
                
                heatmap_fixed = f_freq(std_freqs)

                def find_spectral_nulls(freqs, ts_spectrum):
                    valleys, _ = find_peaks(-ts_spectrum, prominence=10)
                    return freqs[valleys].tolist()
                nulls_found = []
                for i in range(len(target)):
                    nulls_found.append(find_spectral_nulls(proc.f_m/1000, all_ts_specs[i]))
                #print(nulls_found)

                fig, axes = plt.subplots(2, 2, figsize=(16, 12))
                ax_ts = axes[0, 0]    # Haut-Gauche : Spectres TSf
                ax_hm = axes[1, 0]    # Bas-Gauche : Heatmap Angle/Freq
                ax_top = axes[0, 1]   # Haut-Droite : Vue de dessus GPS
                ax_side = axes[1, 1]  # Bas-Droite : Vue de côté (Profondeur)

                for i, t in enumerate(target):
                    color = colors_list[i % 10]
                    if i==0 :
                        ax_ts.plot(proc.f_m/1000, t['ts'], color=color, alpha=0.5, 
                               label=f"P{t['ping']}")
                    else : 
                        ax_ts.plot(proc.f_m/1000, t['ts'], color=color, alpha=0.5, 
                               label=f"P{t['ping']} ({incident_angles[i]:.1f}°)")

                ax_ts.plot(proc.f_m/1000, mean_ts, color='black', lw=3, label='Moyenne 80-100°')
                #ax_ts.set_title(f"Target {idx+1} : Signatures (80-100°)", fontweight='bold')
                ax_ts.set_ylabel("TS (dB)",fontsize=24)
                ax_ts.set_xlabel("Fréquence (kHz)",fontsize=24)
                ax_ts.legend(fontsize=20, ncol=2)
                ax_ts.grid(True, alpha=0.2)

                
                sns.heatmap(heatmap_fixed, ax=ax_hm, cmap='viridis', vmin=-80, vmax=-30,
                            cbar_kws={'label': 'TS (dB)'})
                
                #ax_hm.set_title("Evolution TS : f(Angle, Fréq)", fontweight='bold')
                x_ticks = np.linspace(0, n_freq_std - 1, 5, dtype=int)
                ax_hm.set_xticks(x_ticks + 0.5)
                ax_hm.set_xticklabels(np.round(std_freqs[x_ticks], 1))
                y_ticks = np.linspace(0, n_angle_std - 1, 10, dtype=int)
                ax_hm.set_yticks(y_ticks + 0.5)
                ax_hm.set_yticklabels(np.round(std_angles[y_ticks], 1))
                
                ax_hm.set_ylabel("Angle d'incidence (°)",fontsize=24)
                ax_hm.set_xlabel("Fréquence (kHz)",fontsize=24)


                for i in range(len(xs)):
                    color = colors_list[i % 10]
                    ax_top.plot(xb_target[i], yb_target[i], 's', color=color, ms=5, alpha=0.6)
                    ax_top.plot([xb_target[i], xs[i]], [yb_target[i], ys[i]], ':', color='gray', alpha=0.3)
                    ax_top.plot(xs[i], ys[i], 'o', color=color, ms=10, mec='k', zorder=5)
                #ax_top.set_title("Position GPS (80-100°)", fontweight='bold')
                ax_top.set_aspect('equal', 'datalim')
                ax_top.grid(True, alpha=0.2)


                ax_side.plot(dist_fish_2d, zs, '-', color='blue', alpha=0.2, zorder=1)

                for i in range(len(xs)):
                    ax_side.plot(dist_fish_2d[i], zs[i], 'o', color=colors_list[i % 10], ms=10, mec='k', zorder=5)
                ax_side.invert_yaxis()
                #ax_side.set_title("Profil de Profondeur (Vue de côté)", fontweight='bold')
                ax_side.set_xlabel("Distance Horizontale 2D (m)",fontsize=24)
                ax_side.set_ylabel("Profondeur (m)",fontsize=24)
                ax_side.grid(True, alpha=0.2)

                plt.tight_layout()
                
                
                if not is_rejected :
                    plt.savefig(save_path, dpi=150)
                    plt.close()
                    all_targets_features.append({
                        'depth_m': np.mean([d['z'] for d in target]),
                        'ts_spectra': all_ts_specs,              # Spectres complets pour chaque ping
                        'ts_peak_max': ts_at_peak_db,      # TS au pic de la vessie
                        'tilt_angle': tilt_angle_peak,    # Angle au pic
                        'velocity_m_s': velocity_mps,
                        'ping_count': len(target),
                        'heatmap_flat': heatmap_fixed.flatten()
                    })

                    plot_spectra(target, proc.f_m, os.path.join(file_output_dir, f"Target_{idx+1:03d}_spectra.png"), incident_angles)

        xb = proc.x_boat_gps
        yb = proc.y_boat_gps
        dist_inter_pings = np.sqrt(np.diff(xb)**2 + np.diff(yb)**2)
        profondeurs_fond = r_n[indices_fond]
        angle_ouverture_deg = 7.0 
        angle_ouverture_rad = np.radians(angle_ouverture_deg)
        W_fauchee = 2 * profondeurs_fond * np.tan(angle_ouverture_rad / 2)
        
        S_fauchee = W_fauchee * profondeurs_fond / 2
        print(dist_inter_pings.shape, S_fauchee.shape)
        S_moyennes = (S_fauchee[:-1] + S_fauchee[1:]) / 2
        volume_fichier_m3 = np.sum(dist_inter_pings * S_moyennes)
        
        # 5. Calcul de la densité (nb / 1000 m3)
        nb_poissons = len(all_targets_features)
        densite_1000m3 = (nb_poissons / volume_fichier_m3) * 1000 if volume_fichier_m3 > 0 else 0
        
        # Stockage pour la moyenne finale (à initialiser en haut du script : resultats_collecte = [])
        resultats_globaux.append({

            'fichier': file_id,
            'nb_poissons': nb_poissons,
            'volume': volume_fichier_m3,
            'densite': densite_1000m3
        })
        
        print(f"Volume : {volume_fichier_m3:.1f} m3 | Poissons : {nb_poissons} | Densité : {densite_1000m3:.3f}")

        if all_targets_features:
            df_temp = pd.DataFrame(all_targets_features)
            save_file = os.path.join(output_base_dir, "master_features.pkl")
            
            if os.path.exists(save_file):
                # On charge l'existant, on ajoute le nouveau, on sauvegarde
                df_old = pd.read_pickle(save_file)
                df_temp = pd.concat([df_old, df_temp], ignore_index=True)
            
            df_temp.to_pickle(save_file)
            
            all_targets_features.clear()

    finally:
        if 'proc' in locals(): del proc
        gc.collect()

for f in raw_files:
    try: process_file(f)
    except Exception as e: print(f"Erreur {f}: {e}")

# À LA FIN DU SCRIPT (après la boucle for)
df_final = pd.DataFrame(resultats_globaux)
df_final.to_csv(os.path.join(output_base_dir, "densites_par_fichier.csv"), index=False)

moyenne_densite = df_final['densite'].mean()
densite_ponderee = (df_final['nb_poissons'].sum() / df_final['volume'].sum()) * 1000

print("\n" + "="*30)
print(f"SYNTHÈSE TOTALE ({date.upper()})")
print(f"Nombre total de fichiers : {len(df_final)}")
print(f"Densité moyenne par fichier : {moyenne_densite:.4f} ind/1000m3")
print(f"Densité réelle pondérée : {densite_ponderee:.4f} ind/1000m3")
print("="*30)
