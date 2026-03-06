import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
import os
from scipy.stats import skew

# --- 1. CONFIGURATION ET CHARGEMENT ---
date = "octobre" # ou "octobre"



freq = "70kHz"
output_dir = f"Analyse_Massive_Cibles_{freq}_{date}"
data_path = f"{output_dir}/master_features.pkl"
os.makedirs(output_dir, exist_ok=True)

df = pd.read_pickle(data_path)

# Fonction pour le TS Médian de filtrage
def get_median_ts(specs):
    if not isinstance(specs, np.ndarray): return np.nan
    mean_linear = np.nanmean(10**(specs / 10.0), axis=0)
    mean_db = 10 * np.log10(mean_linear)
    return np.nanmedian(mean_db)

df['ts_median_target'] = df['ts_spectra'].apply(get_median_ts)

# FILTRAGE : On garde les cibles crédibles (-60 à -30 dB)
df = df[(df['ts_median_target'] > -60) & (df['ts_median_target'] < -30)].copy()
print(f"Cibles après filtrage : {len(df)}")

# --- 2. EXTRACTION DES FEATURES MORPHOLOGIQUES ---
def extract_morpho_features(specs_matrix):
    """ Calcule la signature sur l'ensemble de la trace d'un poisson """
    if specs_matrix.ndim == 1: specs_matrix = specs_matrix.reshape(1, -1)
    
    # 1. Moyenne du spectre sur la durée du passage
    mean_spec = np.nanmean(specs_matrix, axis=0)
    
    # 2. Skewness temporel (forme du passage dans le faisceau)
    # Très discriminant pour la morphologie de la vessie
    ts_envelope = 10 * np.log10(np.nanmean(10**(specs_matrix / 10.0), axis=1))
    angular_skew = skew(ts_envelope) if len(ts_envelope) > 3 else 0
    
    # 3. Rugosité fréquentielle
    rugosity = np.mean(np.abs(np.diff(mean_spec)))
    
    # 4. Nombre de pics/creux (interférences)
    valleys, _ = find_peaks(-mean_spec, prominence=1.5)
    
    return [angular_skew, rugosity, len(valleys)]

features_list = []
target_indices = []

for idx, row in df.iterrows():
    specs = row['ts_spectra'] # On attend une matrice [pings x freqs]
    
    if not isinstance(specs, np.ndarray) or specs.ndim < 2:
        continue
    
    # On extrait les features sur la CIBLE entière
    m_feats = extract_morpho_features(specs)
    
    # On ajoute une mesure de stabilité (cohésion entre pings)
    stability = np.mean(np.std(specs, axis=0))
    
    features_list.append(np.append(m_feats, stability))
    target_indices.append(idx)

# Conversion en array pour le GMM
X = np.array(features_list)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

gmm = GaussianMixture(n_components=2, covariance_type='full', n_init=10, random_state=42)
df.loc[target_indices, 'final_cluster'] = gmm.fit_predict(X_scaled)

# --- 4. CALCUL DE LA TAILLE (MOYENNE 45-70 kHz) ---
def calculate_length(row):
    specs = row['ts_spectra']
    cluster = row['final_cluster']
    if not isinstance(specs, np.ndarray) or np.isnan(cluster): return np.nan

    # 1. Moyenne temporelle puis isolation 1ère moitié (45-70 kHz)
    mean_spec = np.nanmean(specs, axis=0) if specs.ndim > 1 else specs
    first_half = mean_spec[len(mean_spec)//3:2*len(mean_spec)//3]
    
    # 2. Moyenne énergétique robuste
    ts_ref = 10 * np.log10(np.nanmean(10**(first_half / 10.0)))
    f_center = 57.5 # Fréquence centrale de la 1ère moitié
    
    # 3. Formules avec dépendance en fréquence
    if cluster == 0: # Cyprinidés
        return 10**((ts_ref + 0.9 * np.log10(f_center) + 62) / 19.1)
    else: # Percidés (Foote)
        return 10**((ts_ref + 0.9 * np.log10(f_center) + 62) / 19.1)

df['length_cm'] = df.apply(calculate_length, axis=1)

# --- 5. VISUALISATION DES SIGNATURES (ANGLE VS FREQ) ---
def plot_signatures(df, output_dir):
    actual_n_freqs = df['ts_spectra'].dropna().iloc[0].shape[-1]
    angles_std = np.linspace(70, 100, 30)
    
    clusters = sorted(df['final_cluster'].dropna().unique())
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    titles = {0: 'Signature : Cyprinidés (C0)', 1: 'Signature : Percidés (C1)'}

    for i, cluster_id in enumerate(clusters):
        subset = df[df['final_cluster'] == cluster_id]
        all_hm_linear = []

        for specs in subset['ts_spectra']:
            if not isinstance(specs, np.ndarray) or specs.ndim < 2: continue
            
            # Interpolation temporelle/angulaire sur 30 points
            specs_linear = 10**(specs / 10.0)
            f_interp = interp1d(np.linspace(70, 100, specs.shape[0]), specs_linear, axis=0, kind='linear', fill_value="extrapolate")
            all_hm_linear.append(f_interp(angles_std))

        if all_hm_linear:
            avg_db = 10 * np.log10(np.nanmean(all_hm_linear, axis=0))
            im = axes[i].imshow(avg_db, aspect='auto', cmap='viridis', vmin=-55, vmax=-45, 
                               origin='lower', extent=[45, 95, 70, 100], interpolation='none')
            #axes[i].set_title(f"{titles[cluster_id]}\n(n={len(all_hm_linear)})", fontweight='bold')
            axes[i].set_xlabel("Fréquence (kHz)", fontsize=24, labelpad=10); axes[0].set_ylabel("Angle d'incidence (°)", fontsize=24, labelpad=10)

            axes[i].tick_params(axis='both', which='major', labelsize=20)
    cbar = fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    cbar.set_label('Target Strength (dB re 1m²)', fontsize=24, labelpad=15)
    cbar.ax.tick_params(labelsize=20)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/Signatures_Angle_Freq.png", dpi=300)

plot_signatures(df, output_dir)

# --- 6. RÉSUMÉ STATISTIQUE ---
print(f"\n--- SYNTHÈSE {date.upper()} ---")
print(df.groupby('final_cluster')[['length_cm', 'depth_m', 'ts_median_target']].mean().round(2))
df.to_pickle(f"{output_dir}/master_features_final.pkl")

# --- 7. TRACÉ DES DISTRIBUTIONS DE TAILLE ---

def plot_size_distributions(df, output_dir, date_str):
    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")
    
    colors = ['#3498db', '#2ecc71'] # Bleu (C0), Vert (C1)
    clusters = sorted(df['final_cluster'].dropna().unique())
    
    for i, cluster_id in enumerate(clusters):
        subset = df[df['final_cluster'] == cluster_id]
        label = "Cyprinidés (C0)" if cluster_id == 0 else "Percidés (C1)"
        
        if not subset['length_cm'].dropna().empty:
            # Tracé du KDE (densité)
            sns.histplot(
                data=subset, x='length_cm', 
                fill=True, alpha=0.6, 
                color=colors[i], label=f"{label} (n={len(subset)})",
                lw=2, binwidth=1

            )
            # Ajout d'une ligne verticale pour la moyenne
            plt.axvline(subset['length_cm'].mean(), color=colors[i], linestyle='--', lw=1.5)

    plt.xlabel("Longueur Totale (cm)", fontsize=24)

    plt.ylabel("Nombre d'individus", fontsize=24)
    plt.legend(fontsize=20)
    
    # Limiter l'axe X pour éviter les queues de distribution aberrantes
    plt.xlim(0, 45)

    plt.tick_params(axis='both', which='major', labelsize=20)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/Distribution_Tailles_{date_str}.png", dpi=300)
    plt.show()

# Appel de la fonction
plot_size_distributions(df, output_dir, date)


# --- 9. GÉNÉRATION DES HISTOGRAMMES DE SYNTHÈSE ---

def plot_final_histograms(df, output_dir, date_str):
    sns.set_style("whitegrid")
    # Création d'une figure à 3 colonnes
    fig, axes = plt.subplots(1, 2, figsize=(22, 7))
    colors = ['#3498db', '#2ecc71'] # C0: Bleu, C1: Vert

    # 2. HISTOGRAMME DES PROFONDEURS (Depth)
    # Inversion de l'axe pour simuler la colonne d'eau
    sns.histplot(data=df, y='depth_m', hue='final_cluster', palette=colors, 
                 kde=False, element="step", ax=axes[0],binwidth=1)
    #axes[0].set_title("Répartition Verticale (Profondeur)", fontsize=14, fontweight='bold')
    axes[0].set_ylabel("Profondeur (m)",fontsize=24)
    axes[0].set_xlabel("Nombre de cibles",fontsize=24)
    axes[0].tick_params(axis='both', which='major', labelsize=20)
    axes[0].invert_yaxis() 

    # 3. HISTOGRAMME DES VITESSES (Speed)
    if 'velocity_m_s' in df.columns:
        sns.histplot(data=df, x='velocity_m_s', hue='final_cluster', palette=colors, 
                     kde=False, element="step", ax=axes[1],binwidth=0.1)

        #axes[1].set_title("Distribution des Vitesses de nage", fontsize=14, fontweight='bold')
        axes[1].set_xlabel("Vitesse (m/s)",fontsize=24)
        axes[1].set_ylabel("Nombre de cibles",fontsize=24)
        axes[1].tick_params(axis='both', which='major', labelsize=20)
    else:
        axes[1].text(0.5, 0.5, "Données de vitesse\nnon disponibles", 
                     ha='center', va='center', fontsize=12, color='red')

    # Ajustement des légendes
    for ax in axes:
        if ax.get_legend():
            ax.get_legend().set_title("Clusters GMM")
            new_labels = ['Cyprinidés (C0)', 'Percidés (C1)']
            for t, l in zip(ax.get_legend().get_texts(), new_labels): t.set_text(l)
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/Histogrammes_Synthese_{date_str}.png", dpi=300)
    plt.show()

# Exécution de la fonction
plot_final_histograms(df, output_dir, date)

# --- 10. CORRÉLATION VITESSE VS TAILLE ---

def plot_velocity_size_correlation(df, output_dir, date_str):
    if 'velocity_m_s' not in df.columns or 'length_cm' not in df.columns:
        print("Données insuffisantes pour le plot de corrélation.")
        return

    plt.figure(figsize=(10, 7))
    sns.set_style("whitegrid")
    
    colors = ['#3498db', '#2ecc71'] # Bleu (C0), Vert (C1)
    
    # Création du scatter plot avec régression linéaire
    for i, cluster_id in enumerate(sorted(df['final_cluster'].dropna().unique())):
        subset = df[df['final_cluster'] == cluster_id]
        label = "Cyprinidés (C0)" if cluster_id == 0 else "Percidés (C1)"
        
        # Plot des points et de la régression
        sns.regplot(
            data=subset, x='length_cm', y='velocity_m_s',
            label=label, color=colors[i],
            scatter_kws={'alpha': 0.5, 's': 30},
            line_kws={'lw': 2}
        )
        
        # Calcul du coefficient de corrélation de Pearson
        corr = subset[['length_cm', 'velocity_m_s']].corr().iloc[0, 1]
        print(f"Corrélation {label} : R = {corr:.2f}")

    plt.title(f"Relation Taille - Vitesse de nage ({date_str.upper()})", fontsize=15, fontweight='bold')
    plt.xlabel("Longueur Totale (cm)", fontsize=12)
    plt.ylabel("Vitesse (m/s)", fontsize=12)
    
    # Optionnel : ajout d'une grille plus fine
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend(title="Classification GMM")
    
    # Ajustement des limites pour la visibilité
    plt.xlim(0, df['length_cm'].max() * 1.1)
    plt.ylim(0, df['velocity_m_s'].max() * 1.1)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/Correlation_Vitesse_Taille_{date_str}.png", dpi=300)
    plt.show()

# Appeler la fonction
plot_velocity_size_correlation(df, output_dir, date)

import seaborn as sns
import matplotlib.pyplot as plt

def plot_joint_distribution(df, output_dir, date_str):
    # Filtrage des NaNs pour éviter les erreurs de tracé
    plot_df = df.dropna(subset=['length_cm', 'velocity_m_s', 'final_cluster'])
    
    # Définition de la palette (C0: Bleu, C1: Vert)
    palette = {0: '#3498db', 1: '#2ecc71'}
    labels = {0: 'Cyprinidés (C0)', 1: 'Percidés (C1)'}

    # Création du JointGrid
    g = sns.JointGrid(data=plot_df, x="length_cm", y="velocity_m_s", hue="final_cluster", palette=palette)

    # 1. Tracé central : Nuage de points + Courbes de densité (KDE 2D)
    g.plot_joint(sns.scatterplot, s=40, alpha=0.4, edgecolor='w')
    g.plot_joint(sns.kdeplot, alpha=0.5, levels=4, linewidths=1.5)


    # 2. Tracés marginaux (Hauts et Droite) : Histogrammes empilés
    g.plot_marginals(sns.histplot, kde=False, alpha=0.3, element="step")

    # Personnalisation des axes et titres
    g.ax_joint.set_xlabel("Longueur Totale (cm)", fontsize=12, fontweight='bold')
    g.ax_joint.set_ylabel("Vitesse de nage (m/s)", fontsize=12, fontweight='bold')
    
    # Légende personnalisée
    handles, _ = g.ax_joint.get_legend_handles_labels()
    g.ax_joint.legend(handles, [labels[0], labels[1]], title="Groupes GMM", loc='upper left')

    plt.suptitle(f"Analyse Croisée Taille vs Vitesse - {date_str.upper()}", y=1.02, fontsize=16, fontweight='bold')
    
    # Sauvegarde
    plt.savefig(f"{output_dir}/JointPlot_Taille_Vitesse_{date_str}.png", dpi=300, bbox_inches='tight')
    plt.show()

# Exécution
plot_joint_distribution(df, output_dir, date)