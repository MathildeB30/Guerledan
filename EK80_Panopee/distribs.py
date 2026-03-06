import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- 1. CHARGEMENT ---
pkl_path = "Analyse_Massive_Cibles_70kHz_fevrier/master_features.pkl"
df = pd.read_pickle(pkl_path)

# --- 2. FONCTION DE MOYENNAGE SUR LE PREMIER TIERS ---
def mean_first_third(spectrum):
    """
    Prend un array de TS (spectre), calcule le premier tiers
    et en extrait la moyenne énergétique (linéaire puis retour en dB).
    """
    if isinstance(spectrum, (list, np.ndarray)):
        spectrum = np.array(spectrum)
        # On définit l'index de coupure (1/3 du spectre)
        cut_off = len(spectrum) // 3
        first_third = spectrum[cut_off:2*cut_off]
        
        if len(first_third) > 0:
            # Pour une moyenne acoustique correcte, on passe en linéaire, on moyenne, puis retour en dB
            # Formule : 10 * log10( moyenne( 10^(TS/10) ) )
            lin_val = 10**(first_third / 10)
            mean_lin = np.nanmean(lin_val)
            return 10 * np.log10(mean_lin)
    
    return np.nan # Si la donnée n'est pas un spectre

# --- 3. TRAITEMENT DES DONNÉES ---
# Supposons que votre colonne de spectre s'appelle 'ts_spectrum' ou 'ts_mean' (si elle contient des arrays)
col_source = 'ts_spectra' if 'ts_spectra' in df.columns else 'ts_mean'


print(f"Calcul du TS moyen sur le premier tiers de la colonne : {col_source}")
df['ts_mean_1_3'] = df[col_source].apply(mean_first_third)

# Suppression des échos vides ou erreurs de calcul
df = df.dropna(subset=['ts_mean_1_3'])

# --- 4. CALCUL DE LA TAILLE (LOVE 70kHz) ---
# L = 10^((TS + 68) / 20)
df['length_cm'] = 10**((df['ts_mean_1_3'] + 62 + 0.9*np.log10(57)) / 19.1)

# --- 5. VISUALISATION ---
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
sns.set_style("whitegrid")

# Distribution du TS (1/3 spectre)
sns.histplot(df['ts_mean_1_3'], kde=True, ax=axes[0], color='blue')
axes[0].set_title('Distribution TS')
axes[0].set_xlabel('TS Moyen (dB)')

# Distribution des Tailles
sns.histplot(df['length_cm'], kde=True, ax=axes[1], color='green')
axes[1].set_title('Tailles Estimées à Guerlédan (cm)')
axes[1].set_xlabel('Longueur (cm)')

plt.tight_layout()
plt.savefig("Analyse_Massive_Cibles_70kHz_fevrier/distribs.png")




print(f"Analyse terminée sur {len(df)} cibles.")
print(df[['ts_mean_1_3', 'length_cm']].describe().round(2))

# --- 6. IDENTIFICATION DES GROS INDIVIDUS ---

# 1. On trie le DataFrame pour voir les 20 plus grands
print("\n" + "="*50)
print("LES 20 PLUS GRANDS INDIVIDUS DÉTECTÉS")
print("="*50)
# On sélectionne les colonnes utiles (ajoute 'final_cluster' si tu as fait le clustering)
cols_to_show = ['ts_mean_1_3', 'length_cm']
if 'depth_m' in df.columns: cols_to_show.append('depth_m')
if 'final_cluster' in df.columns: cols_to_show.append('final_cluster')

top_20 = df.sort_values(by='length_cm', ascending=False).head(20)
print(top_20[cols_to_show].round(2))

# 2. Filtrage par seuil (ex: individus de plus de 40 cm)
seuil_taille = 20
big_fish = df[df['length_cm'] > seuil_taille]

print("\n" + "="*50)
print(f"STATISTIQUES SUR LES INDIVIDUS > {seuil_taille} cm")
print("="*50)
print(f"Nombre total trouvé : {len(big_fish)}")

if not big_fish.empty:
    print(f"Taille moyenne : {big_fish['length_cm'].mean():.2f} cm")
    if 'depth_m' in big_fish.columns:
        print(f"Profondeur moyenne : {big_fish['depth_m'].mean():.2f} m")
    
    # Affichage des 5 plus gros parmi ceux-là
    print("\nTop 5 des plus gros spécimens :")
    print(big_fish[cols_to_show].sort_values(by='length_cm', ascending=False).head(20).round(2))
else:
    print(f"Aucun poisson n'a été détecté au-dessus de {seuil_taille} cm.")