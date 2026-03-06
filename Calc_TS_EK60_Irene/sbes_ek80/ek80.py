""" 
Décodage d'un fichier raw d'un EK80

LIMITATION : pour l'instant on suppose qu'il n'y a qu'un capteur

TODO : assert monovoie
"""

import os
import math
import struct
import datetime
import sys
import numpy as np

from . import sbes_Kongsberg
from . import seawater_sound_absorption
from .NMEA import NMEA_decode

class EK80File:
    def __init__(self, directory_in, filename, q_debug=None):
        """
        Lecture d'un ou plusieurs fichiers de sondeurs d'EK80
        """
        # sbes_Kongsberg ne permet que des listes et des tuples
        if type(filename not in (list, tuple)):
            filename = (filename,)

        # q_debug should be a list or a tuple
        if type(q_debug not in (list, tuple)):
            q_debug = (q_debug,)
        
        self.q_debug = q_debug

        # Nombre de ping lus dans le fichier
        self.n0_ping = 0
        
        # Ouverture du fichier
        self.ek80 = sbes_Kongsberg.EK80(directory_in, filename, None,
                                        self.q_debug)

        # First file address since last ping
        self.last_ping_address = 0
        
    def read_file(self, ping0, n_ping_, q_angle, q_normalize):
        """ 
        Lecture d'une partie du fichier filename 
        
        ping0: int : numéro du premier ping à lire
                     ping0=None : on lit à partir du ping 0
        n_ping_: int : nombre de pings à lire 
                      n_ping=None : lire tous les pings
        q_angle: bool : si True : on lit les angles d'arrivées
        q_normalize: si "TS" : on convertit les niveaux en TS
                     si "SV" : on convertit les niveaux en SV
                     si "None" : on convertit en puissance reçue
        """

        q_eof = False                    
        if n_ping_ is None:
            # On suppose que l'on atteindra jamais la valeur suivante...
            ping1 = 2 ** 31
        else:
            ping1 = ping0 + n_ping_

        
        # Initialisation des différentes structures
        n0_ping = 1024
        n_ping = 0
        n0_range = 0
        self.t_ping = np.empty((n0_ping,), np.float64)
        self.amplitude_ping = np.empty((n0_ping, n0_range), np.float32)
        if q_angle == True:
            self.angle_ping = np.empty((n0_ping, n0_range, 2), np.float32)
                
        n0_pos = 1024
        n_pos = 0
        self.t_pos = np.empty((n0_pos,), np.float64)
        self.pos_pos = np.empty((n0_pos, 3), np.float64)

        n0_course = 1024
        n_course = 0
        self.t_course = np.empty((n0_course,), np.float64)
        self.speed_course = np.empty((n0_course,), np.float32)
        self.heading_course = np.empty((n0_course,), np.float32)
    
        n0_att = 1024
        n_att = 0
        self.t_att = np.empty((n0_att,), np.float64)
        self.att_att = np.empty((n0_att, 3), np.float32)
        
        n0_heave = 1024
        n_heave = 0
        self.t_heave = np.empty((n0_heave,), np.float64)
        self.heave_heave = np.empty((n0_heave,), np.float32)

            
        # On se place au premier ping : lecture à blanc des pings précédents
        # On garde quand même les variables xml
        #--------------------------------------------------------------------
        if ping0 is None:
            ping0 = 0

        if ping0 - self.n0_ping < 0:
            # On repart de 0
            print("File rewinded to continue process")
            self.ek80.rewind()
            self.n0_ping = 0

        while self.n0_ping < ping0:
            try:
                type_ = self.ek80.read_packet()
                if type_ == b"XML0":
                    # On ne conserve que les dernières valeurs lues
                    self.process_xml(self.ek80.xml_type)

                elif type_ == b"RAW3":
                    self.n0_ping += 1
                    self.last_ping_address = self.ek80.f.tell()
                    
            except EOFError:
                print("Not enough ping to start at {:d}th ping (file {:s})"\
                      .format(ping0, self.ek80.filelist[0]))
                q_eof = True
                break
            
        # Lecture des pings que l'on conserve
        #------------------------------------

        # On repart de la dernière adresse lue
        self.ek80.f.seek(self.last_ping_address, os.SEEK_SET)
        n_ping = 0

        
        while 1:
            try:
                type_ = self.ek80.read_packet()
                
                if type_ == b"XML0":
                    self.process_xml(self.ek80.xml_type)
        
                elif type_ == b"MRU0":
                
                    if n_att >= n0_att:
                        n0_att *= 2
                        self.t_att = np.resize(self.t_att, (n0_att,))
                        self.att_att = np.resize(self.att_att, (n0_att, 3))
                        
                    if n_heave >= n0_heave:
                        n0_heave *= 2
                        self.t_heave = np.resize(self.t_heave, (n0_heave,))
                        self.heave_heave = np.resize\
                            (self.heave_heave, (n0_heave,))

                    self.t_att[n_att] = self.ek80.time
                    self.att_att[n_att,:] = np.array\
                        ((self.ek80.mru["heading"], self.ek80.mru["pitch"],
                          self.ek80.mru["roll"]))

                    n_att += 1

                    self.t_heave[n_heave] = self.ek80.time
                    self.heave_heave[n_heave] = self.ek80.mru["heave"]
                    n_heave += 1

                elif type_ == b"RAW3":

                    # Continue t'on la lecture ?
                    if self.n0_ping >= ping1:
                        # Stop recurrence
                        break
                    
                    # Mise à jour de l'allocation
                    if n_ping >= n0_ping:
                        n0_ping *= 2
                        self.t_ping = np.resize(self.t_ping, (n0_ping,))
                        self.amplitude_ping = np.resize\
                            (self.amplitude_ping, (n0_ping, n0_range))
                        if q_angle == True:
                            self.angle_ping = np.resize\
                                (self.angle_ping, (n0_ping, n0_range, 2))

                    # On met à jour aussi le range si nécessaire
                    n_range= self.ek80.raw["power"].size      
                
                    if n_range > n0_range:
                        dn = n_range - n0_range
                        n0_range = n_range
                        self.amplitude_ping = np.concatenate\
                            ((self.amplitude_ping,
                              np.nan * np.ones((n0_ping, dn), np.float32)), 1)
                        if q_angle == True:
                            self.angle_ping = np.concatenate\
                                ((self.angle_ping,
                                  np.nan * np.ones((n0_ping, dn, 2),
                                                   np.float32)), 1)
                            
                    # Chargement des données
                    self.t_ping[n_ping] = self.ek80.time
                    self.amplitude_ping[n_ping, :n_range]\
                        = self.ek80.raw["power"]
                    self.amplitude_ping[n_ping, n_range:] = np.nan
                    if q_angle == True:
                        self.angle_ping[n_ping, :n_range, :]\
                            = self.ek80.raw["phase"]
                        self.angle_ping[n_ping, n_range:, :] = np.nan

                    # Vérification de la non-modification de paramètres critique
                    #power = float(self.ek80.parameter["Channel"]["TransmitPower"])
                    #frequency = float(self.ek80.parameter["Channel"]["Frequency"])
                    #pulse_duration\
                    #    = float(self.ek80.parameter["Channel"]["PulseDuration"])
                    #sample_interval\
                    #    = float(self.ek80.parameter["Channel"]["SampleInterval"])

                    #try:
                    #    sound_velocity\
                    #    = float(self.ek80.parameter["Channel"]["SoundVelocity"])
                    #except KeyError:
                    #    sound_velocity = self.average_sound_velocity
                    
                    if n_ping == 0:
                        # Initialisation
                        self.frequency_ = self.frequency
                        self.power_ = self.power
                        self.pulse_duration_ = self.pulse_duration
                        self.sample_interval_ = self.sample_interval
                        self.sound_velocity_ = self.sound_velocity
                        
                    else:
                        # Test des vérifications
                        if self.frequency_ != self.frequency:
                            print("frequency has changed (from {:f} to {:f}"\
                                  .format(self.frequency_, self.frequency))
                            break
                        
                        if self.power_ != self.power:
                            print("power has changed (from {:f} to {:f}"\
                                  .format(self.power_, self.power))
                            break

                        #print(self.pulse_duration_, self.pulse_duration)
                        if self.pulse_duration_ != self.pulse_duration:
                            print("pulse duration has changed"\
                                  +" (from {:f} to {:f}"\
                                  .format(self.pulse_duration_,\
                                          self.pulse_duration))
                            break
                        if self.sample_interval_ != self.sample_interval:
                            print("sample interval has changed"\
                                  +" (from {:f} to {:f}"\
                                  .format(self.sample_interval_,
                                          self.sample_interval))
                            break
                        if self.sound_velocity_ != self.sound_velocity:
                            print("sound velocity has changed"\
                                  + " (from {:f} to {:f}"\
                                  .format(self.sound_velocity_,
                                          self.sound_velocity))
                            break

                    self.last_ping_address = self.ek80.f.tell()
                    n_ping += 1
                    self.n0_ping += 1

                
                elif type_ == b"NME0":

                    #print(self.ek80.nmea)
                    #import pdb;pdb.set_trace()
                    
                    frame = NMEA_decode(self.ek80.nmea)
                    if frame["type"] == "GGA":
                        
                        if n_pos >= n0_pos:
                            n0_pos *= 2
                            self.t_pos = np.resize(self.t_pos, (n0_pos,))
                            self.pos_pos = np.resize(self.pos_pos, (n0_pos, 3))

                        self.t_pos[n_pos] = self.ek80.time
                        self.pos_pos[n_pos, 0] = frame["longitude"] * 180. / np.pi
                        self.pos_pos[n_pos, 1] = frame["latitude"] * 180. / np.pi
                        self.pos_pos[n_pos, 2] = frame["altitude"]
                        n_pos += 1
                    
                    elif frame["type"] == "VTG":
                        if n_course >= n0_course:
                            n0_course *= 2
                            self.t_course = np.resize(self.t_course, (n0_course,))
                            self.speed_course = np.resize(self.speed_course, (n0_course,))
                            self.heading_course = np.resize(self.heading_course, (n0_course,))
                        self.t_course[n_course] = self.ek80.time
                        self.speed_course[n_course] = frame["speed"]
                        self.heading_course[n_course] = frame["true heading"]\
                            * np.pi / 180.
                        
                        n_course += 1
                    
            except EOFError:
                # Fin au moins provisoire de la lecture du fichier
                q_eof = True
                break

        # On clôt la lecture
        if n_ping > 0:
            # Fin du fichier (ou erreur détectée)
            self.t_ping = np.resize(self.t_ping, (n_ping,))
            self.amplitude_ping = np.resize(self.amplitude_ping,
                                            (n_ping, n0_range))
            if q_angle == True:
                self.angle_ping = np.resize(self.angle_ping,
                                            (n_ping, n0_range, 2))
            else:
                self.angle_ping = None

            if n_att > 0:
                # attitude
                self.t_att = np.resize(self.t_att, (n_att,))
                self.att_att = np.resize(self.att_att, (n_att, 3))
            else:
                self.t_att = None
                self.att_att = None

            if n_heave > 0:
                self.t_heave = np.resize(self.t_heave, (n_heave,))
                self.heave_heave = np.resize(self.heave_heave, (n_heave,))
            else:
                self.t_heave = None
                self.heave_heave = None

            if n_pos > 0:
                self.t_pos = np.resize(self.t_pos, (n_pos,))
                self.pos_pos = np.resize(self.pos_pos, (n_pos, 3))
            else:
                self.t_pos = None
                self.pos_pos = None

            if n_course > 0:
                self.t_course = np.resize(self.t_course, (n_course,))
                self.speed_course = np.resize(self.speed_course, (n_course,))
                self.heading_course = np.resize(self.heading_course, (n_course,))
            else:
                self.t_course = None
                self.speed_course = None
                self.heading_course = None            

        # Date de référence pour tous les temps
        self.reference_time = self.ek80.ref_time
        
        # distance max en temps
        self.time_range = n_range * self.sample_interval

        # Fin de la correction des angles
        if q_angle == True:
            self.angle_ping[..., 0] = self.angle_ping[..., 0]\
                / self.angle_sensitivity_alongship
            self.angle_ping[...,1] = self.angle_ping[..., 1]\
                / self.angle_sensitivity_athwartship
            
        path_cal = "EK60_Ulysse/Fevrier/DATA/CalibrationDataFile-D20260205-T175310.xml" 
        
        if os.path.exists(path_cal):
            self.apply_external_calibration(path_cal)

        if q_normalize == "TS":
            self.TS_compute()
        elif q_normalize == "SV":
            self.Sv_compute()
                
        return q_eof
    
    def apply_external_calibration(self, xml_path):
        """ Charge et applique des paramètres de calibration depuis un XML externe """
        import xml.etree.ElementTree as ET
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            # On cherche les valeurs dans le XML (adapter les balises si nécessaire)
            # Note : Simrad utilise souvent des structures complexes, on utilise .find()
            ext_gain = float(root.find(".//Gain").text)
            ext_sa = float(root.find(".//SaCorrection").text)
            
            # On écrase les listes de gain/SaCorr pour toutes les durées d'impulsion
            # Ou on peut cibler une durée précise si besoin.
            self.gain_list = [ext_gain] * len(self.gain_list)
            self.sa_correction_list = [ext_sa] * len(self.sa_correction_list)
            
            print(f"--- Calibration externe appliquée : G={ext_gain}, SaCorr={ext_sa} ---")
        except Exception as e:
            print(f"Erreur lors de la lecture du XML de calibration : {e}")
    
    def process_xml(self, xml_type):
        """ 
        Extraction des paramètres pertinents provenant
            des trames xml
        """
        if xml_type == "Configuration":
            # Paramètres de configuration
            self.equivalent_beam_angle\
                = float(self.ek80.configuration["Transceivers"]\
                        ["Transceiver"]["Channels"]["Channel"]\
                        ["Transducer"]["EquivalentBeamAngle"])
            
            self.angle_sensitivity_alongship\
                = float(self.ek80.configuration["Transceivers"]\
                        ["Transceiver"]["Channels"]["Channel"]\
                        ["Transducer"]["AngleSensitivityAlongship"])
            
            self.angle_sensitivity_athwartship\
                = float(self.ek80.configuration["Transceivers"]\
                        ["Transceiver"]["Channels"]["Channel"]\
                        ["Transducer"]["AngleSensitivityAthwartship"])
            
            self.beamwidth_alongship\
                = float(self.ek80.configuration["Transceivers"]\
                        ["Transceiver"]["Channels"]["Channel"]\
                        ["Transducer"]["BeamWidthAlongship"])
            
            self.beamwidth_athwartship\
                = float(self.ek80.configuration["Transceivers"]\
                        ["Transceiver"]["Channels"]["Channel"]\
                        ["Transducer"]["BeamWidthAthwartship"])
            
            # Tableau des gains et SaCorr par longueur de pulse
            txt = self.ek80.configuration["Transceivers"]\
                ["Transceiver"]["Channels"]["Channel"]\
                ["PulseDuration"]
            
            self.pulse_duration_list\
                = [ float(x) for x in txt.split(";") ]
            
            txt = self.ek80.configuration["Transceivers"]\
                ["Transceiver"]["Channels"]["Channel"]\
                ["Transducer"]["Gain"]
            self.gain_list\
                = [ float(x) for x in txt.split(";") ]

            txt = self.ek80.configuration["Transceivers"]\
                ["Transceiver"]["Channels"]["Channel"]\
                ["Transducer"]["SaCorrection"]
            self.sa_correction_list\
                = [ float(x) for x in txt.split(";") ]

        elif xml_type == "Environment":            
            # Données d'environnement
            self.average_sound_velocity\
                = float(self.ek80.environment["SoundSpeed"])
            
            self.temperature\
                = float(self.ek80.environment["Temperature"])

            self.salinity\
                = float(self.ek80.environment["Salinity"])

            self.pH\
                = float(self.ek80.environment["Acidity"])

            self.depth\
                = float(self.ek80.environment["Depth"])

            # Profil de célérité    
            txt = self.ek80.environment["SoundVelocityProfile"]
            
            svp = np.fromstring(txt, np.float32,-1,";")
            svp.shape = (-1, 2)
            self.sound_velocity_profile = svp

        elif xml_type == "Parameter":
            
            # Vérification de la non-modification de paramètres critiques
            self.power = float(self.ek80.parameter["Channel"]["TransmitPower"])
            self.frequency = float(self.ek80.parameter["Channel"]["Frequency"])
            self.pulse_duration\
                = float(self.ek80.parameter["Channel"]["PulseDuration"])
            self.sample_interval\
                = float(self.ek80.parameter["Channel"]["SampleInterval"])
            try:
                self.sound_velocity\
                    = float(self.ek80.parameter["Channel"]["SoundVelocity"])
            except KeyError:
                self.sound_velocity = self.average_sound_velocity
            
    def process_NMEA(self, nmea):
        frame = NMEA_decode(self.ek80.nmea)
        if frame["type"] == "GGA":

            self.q_pos = True
            self.t_pos_ = self.ek80.time
            self.pos_pos_ \
                = np.array((frame["longitude"] * 180. / np.pi,
                            frame["latitude"] * 180. / np.pi, 
                            frame["altitude"]))

        elif frame["type"] == "VTG":
            self.q_course = True
            self.t_course_ = self.ek80.time
            self.speed_course_ = frame["speed"]
            self.course_course_ = frame["course"]


    def process_MRU0(self, mru):
        self.q_att = True
        self.t_att_ = self.ek80.time
        self.att_att_ = np.array\
            ((self.ek80.mru["heading"], self.ek80.mru["pitch"],
              self.ek80.mru["roll"]))

        self.q_heave = True
        self.t_heave_ = self.ek80.time
        self.heave_heave_ = self.ek80.mru["heave"]
        
    def TS_compute(self):
        """ Utilisation de la formule de echoview pour corriger le TS """
        # Distances
        nr = self.amplitude_ping.shape[1]
        dr = self.sample_interval  * self.average_sound_velocity / 2.
    
        r = np.arange(nr) * dr
        # On modifie r[0] pour éviter log(0) dans TL
        r[0] = dr
        
    
        # Coefficient d'absorption
        alpha = seawater_sound_absorption.Francois_Garrison\
            (self.frequency, self.pH, self.salinity, self.temperature, self.depth)

        # Sélection du gain et du Sa corr
        try:
            idx = self.pulse_duration_list.index\
                (self.pulse_duration)
        except:
            print("Error in gain and sa correction correction:"\
                  + "pulse duration: {:g} not found!".format\
                  (self.pulse_duration))
            sys.exit()
        
        lbd = self.sound_velocity / self.frequency
    
        TL = 40. * np.log10(r) + 2 * alpha *r
        SL = 10. * np.log10(self.power)
        G = self.gain_list[idx]
        DI_tx = G - 10 * np.log10(4. * np.pi)
        DI_rx = G + 10 * np.log10(lbd ** 2 /(4. * np.pi))
        
        self.amplitude_ping += TL - SL - DI_tx - DI_rx

    def Sv_compute(self):
        """ Utilisation de la formule de echoview pour corriger le TS """
        
        # Distances
        nr = self.amplitude_ping.shape[1]
        dr = self.sample_interval * self.average_sound_velocity / 2.
    
        r = np.arange(nr) * dr
        # On modifie r[0] pour éviter log(0) dans TL
        r[0] = dr
        
        # Coefficient d'absorption
        alpha = seawater_sound_absorption.Francois_Garrison\
            (self.frequency, self.pH, self.salinity, self.temperature, self.depth)

        # Sélection du gain et du Sa corr
        try:
            idx = self.pulse_duration_list.index\
                (self.pulse_duration)
        except:
            print("Error in gain and sa correction correction:"\
                  + "pulse duration: {:g} not found!".format\
                  (self.pulse_duration))
            sys.exit()

        c = self.sound_velocity
        lbd = self.sound_velocity / self.frequency
        tau = self.pulse_duration
        
        TL = 40. * np.log10(r) + 2 * alpha *r
        SL = 10. * np.log10(self.power)
        G = self.gain_list[idx]
        Psi = self.equivalent_beam_angle
        Sa_corr = 2. * self.sa_correction_list[idx]
        DI_tx = G - 10 * np.log10(4. * np.pi)
        DI_rx = G + 10 * np.log10(lbd ** 2 /(4. * np.pi))

        V = 10 * np.log10(c * tau / 2) + Psi + 20. * np.log10(r) + 2 * Sa_corr
        
        self.amplitude_ping\
            += (TL - SL - DI_tx - DI_rx - V)[None, :]

    
