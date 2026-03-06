Post-traitement des données sondeur Simrad EK80

Ce dépôt contient un ensemble d'outils Python dédiés à l'extraction, à la normalisation et à la conversion des données acoustiques brutes issues des sondeurs scientifiques Simrad EK60.

Description des composants

m_export_EK80.py : Script principal pilotant le flux de données. Il permet d'exporter les pings acoustiques vers des formats compatibles avec les environnements d'analyse de données standards tels que Matlab (.mat) ou NumPy (.npz).

Module ek80.py (sbes_ek80) : Ce module assure le calcul des indicateurs acoustiques TS (Target Strength) et Sv (Volume Scattering). Il intègre les paramètres de calibration ainsi que les profils de célérité (SVP) pour garantir la précision des données en post-traitement.

Module NMEA.py : Responsable de la gestion des données de navigation. Il extrait les trames NMEA et effectue une interpolation temporelle pour synchroniser la position (GPS), la vitesse, le cap et les mouvements de la plateforme (roulis, tangage, pilonnement) avec chaque émission acoustique.

Configuration et utilisation

Le comportement du programme est défini par le fichier de configuration EK80.yaml. L'utilisateur doit y renseigner les paramètres suivants avant l'exécution :

Chemins d'accès : Définition des répertoires d'entrée (fichiers .raw) et de sortie.

Type de données : Choix du niveau de normalisation souhaité (TS ou Sv).

Paramètres d'export : Sélection du format de fichier final (.mat ou .npz) et définition du radical des noms de fichiers.

Flux de traitement

Le script suit une logique de traitement rigoureuse :

Lecture du flux binaire Simrad.

Synchronisation des données de mouvement par interpolation sur les pings.

Application des gains de calibration.

Visualisation optionnelle pour validation avant l'exportation finale.
