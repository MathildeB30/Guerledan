#!/usr/bin/python3
"""

Export de données d'un fichier EK80

Version 1.0 : MLG  2024/04/19

Execution :
python m_export_data <fichier.raw> 
    exportation de tous les pings du fichier
    
python m_export_data <fichier.raw> <ping0>
    exportation du fichier à partir du <ping0> jusqu'à la fin


python m_export_data <fichier.raw> <ping0> <nping>
    exportation du fichier sur <nping> à partir du ping <ping0>


"""
import sys
import os
import datetime

import argparse
import numpy as np
import matplotlib.pyplot as plt
import yaml

import scipy.io
import scipy.interpolate

import sbes_ek80




if __name__ == "__main__":

    # Définition des valeurs d'entrées
    #-----------------------------------

    # Il y a trois niveaux d'entrée des paramètres modifiables du programme
    # Par ordre de priorité croissante :
    #     les valeurs par défaut écrite en dur ci-dessous
    #     les valeurs écrites dans le fichier yaml : EK80.yaml
    #     les valeurs entrées en paramètres
    
    # Liste des valeurs par défaut
    #------------------------------
    # Normalisation des données en TS, SV ou données brutes (puissance reçue)
    q_normalize = "SV"  # "TS", "SV" ou None

    # Visualisation des données
    q_visu = False  # False ou True

    # Palette de couleur pour la visualisation de l'amplitude
    colormap = "simrad" # Une des couleurs de Matplotlib ou "simrad"
    colorbar_lim = None  # None, (None, v_max), (v_min, None), (v_min, v_max)

    # Format de sortie
    q_export = ".mat"  # format ".npz" ou format ".mat"

    # Affichage d'informations de lecture du fichier
    # None, "xml set", "xml", "xml raw", "nmea", "nmea set"
    q_debug = None
    
    # Répertoire d'entrée par défaut
    directory_in = "."
    
    # Répertoire de sortie par défaut
    directory_out = "./EK60_Ulysse/Fevrier/echo_irene/"

    # Lit on les déphasages ?
    q_angle = False

    # Radical du fichier de sortie
    radix_out = "out"
    
    # Liste des fichiers par défaut
    filelist = (
        "EK60ENSTAB-D20240418-T083428.raw",
        "EK60ENSTAB-D20240418-T085052.raw",
        "EK60ENSTAB-D20240418-T090510.raw",
        "EK60ENSTAB-D20240418-T091807.raw",
        "EK60ENSTAB-D20240418-T091817.raw",
	"EK60ENSTAB-D20240418-T092218.raw",
	"EK60ENSTAB-D20240418-T092757.raw",
	"EK60ENSTAB-D20240418-T093549.raw",
	"EK60ENSTAB-D20240418-T094338.raw",
	"EK60ENSTAB-D20240418-T094849.raw",
	"EK60ENSTAB-D20240418-T095558.raw",
	"EK60ENSTAB-D20240418-T100504.raw",
	"EK60ENSTAB-D20240418-T100751.raw",
	"EK60ENSTAB-D20240418-T101229.raw",
	"EK60ENSTAB-D20240418-T101709.raw",
	"EK60ENSTAB-D20240418-T103352.raw",
	"EK60ENSTAB-D20240418-T103919.raw",
	"EK60ENSTAB-D20240418-T104506.raw",
	"EK60ENSTAB-D20240418-T110400.raw",
	"EK60ENSTAB-D20240418-T110852.raw",
	"EK60ENSTAB-D20240418-T112123.raw")

    # Numéro du premier ping à lire dans le fichier
    ping0 = None
    # Nombre de pings à lire par fichier de sortie
    n_ping = None
    # Un seul fichier de sortie ou tous ?
    q_iterate = True

    # Lecture des fichiers yaml
    #----------------------------
    try:
        with open("Calc_TS_EK60_Irene/EK80.yaml", "rt") as f:

            data = yaml.safe_load(f)



            if 'q_normalize' in data.keys():
                q_normalize = data["q_normalize"]
                
            if 'q_visu' in data.keys():
                q_visu = data["q_visu"]
                
            if 'q_angle' in data.keys():
                q_angle = data["q_angle"]

            if 'colormap' in data.keys():
                colormap = data["colormap"]

            if 'colorbar_lim' in data.keys():
                colorbar_lim = data["colorbar_lim"]

            if 'q_export' in data.keys():
                q_export = data["q_export"]

            if 'q_debug' in data.keys():
                q_debug = data["q_debug"]

            if 'directory_in' in data.keys():
                directory_in = data["directory_in"]

            if 'directory_out' in data.keys():
                directory_out = data["directory_out"]

            if 'radix_out' in data.keys():
                radix_out = data["radix_out"]

            if 'filelist' in data.keys():
                filelist = data["filelist"]

            if 'ping0' in data.keys():
                ping0 = data["ping0"]

            if 'n_ping' in data.keys():
                n_ping = data["n_ping"]

            if 'q_iterate' in data.keys():
                q_iterate = data["q_iterate"]
                
    except yaml.YAMLError as exc:
        if hasattr(exc, "problem_mark"):
            mark = exc.problem_mark
            print("Error parsing yaml file at line {:d}, column {:d}"\
                  .format(mark.line, mark.column+1))
            if exc.context != None:
                print(exc.problem_mark)
                print(exc.problem, " ", exc.context)
            else:
                print(exc.problem_mark)
                print(exc.problem)
            sys.exit()


    # Surcharge par les paramètres entrées
    #--------------------------------------
    parser = argparse.ArgumentParser(
        prog='m_export_EK80.py',
        description=
        """
        
        Program which parse a list of EK80 raw data file and allow to export main data
        to ".npz" numpy file format or ".mat" Matlab file format.
        It allows optionally to:\n
               - Scale watercolumn data to TS level (Target Strength Level) or SV level\n
                                  (Volume scattering level)\n
                                - Export attitude and position data at ping time\n
                                - Visualize data\n
        """)
    parser.add_argument('filename', nargs="?",
                        help="filename to process")
    parser.add_argument('radix', nargs="?",
                        help="radix of the output file")
    parser.add_argument('ping0', type=int, nargs="?",
                        help="number of first ping to read")
    parser.add_argument('n_ping', type=int, nargs="?",
                        help="number of ping to read")
    parser.add_argument('-i', '--iterate', action="store_true",
                        required=False,
                        help="if present, batch mode: "
                        + "sequences are read up to end of file by size of n_ping")
    
    args = parser.parse_args()
    if args.filename is not None:
        filelist = [args.filename]
        directory_in = "."
    if args.radix is not None:
        radix_out = args.radix
    if args.ping0 is not None:
        ping0 = args.ping0
    if args.n_ping is not None:
        n_ping = args.n_ping
    if args.iterate == True:
        q_iterate = True

    # Quelques rares formattages
    if ping0 is None:
        ping0 = 0
        
    if type(filelist) not in (list, tuple):
            filelist = [filelist]
        
    # Numéro de figure si l'affichage est demandé
    nro_fig = 0
    
    for filename in filelist:

        # Numéro du batch si iterate
        nro_batch = 1
        ping0_ = ping0
        print("Process file {:s}".format(filename))
        
        # Ouverture du fichier
        ek80 = sbes_ek80.EK80File(directory_in, filename, q_debug)

        while 1:
            # lecture du fichier
            q_eof = ek80.read_file\
                (ping0_, n_ping, q_angle, q_normalize)

            if q_visu == True:
                # Affichage des images
                if colormap == "simrad":
                    colormap = sbes_ek80.get_Simrad_colorbar()

                
                # Image d'amplitude
                if q_normalize == "TS":
                    title = r"TS (dB ref 1m$^{-2}$)"
                elif q_normalize == "SV":
                    title = "SV (dB ref 1m)"
                else:
                    title = "received power (dB)"
                
                
                ff0, ax0 =sbes_ek80.plot_image\
                    (img=ek80.amplitude_ping.T,
                     x_label="ping #",
                     y_label="time (ms)",
                     title=title,
                     nro_fig=nro_fig,
                     xlim=None,
                     ylim=[ek80.time_range * 1000., 0],
                     zlim=colorbar_lim,
                     q_origin="upper",
                     colormap=colormap,
                     q_grid=False,
                     q_interpolation="linear",
                     aspect="auto",
                     q_colorbar=True)
                
                if q_angle:
                    ff1, ax1 =sbes_ek80.plot_image\
                        (img=ek80.angle_ping[...,0].T,
                         x_label="ping #",
                         y_label="time (ms)",
                         title="athwartship angle (degrees)",
                         nro_fig=nro_fig+1,
                         xlim=None,
                         ylim=[ek80.time_range * 1000., 0],
                         zlim=None,
                         q_origin="upper",
                         colormap="jet",
                         q_grid=False,
                         q_interpolation="linear",
                         aspect="auto",
                         q_colorbar=True)
                
                    ff2, ax2 =sbes_ek80.plot_image\
                        (img=ek80.angle_ping[...,1].T,
                         x_label="ping #",
                         y_label="time (ms)",
                         title="longitudinal angle (degrees)",
                         nro_fig=nro_fig+2,
                         xlim=None,
                         ylim=[ek80.time_range * 1000., 0],
                         zlim=None,
                         q_origin="upper",
                         colormap="jet",
                         q_grid=False,
                         q_interpolation="linear",
                         aspect="auto",
                         q_colorbar=True)

                if ek80.t_att is not None:
                    ff4, ax4 = sbes_ek80.plot_xy\
                        (ek80.t_att,
                         ek80.att_att[:,0] * 180. / np.pi,
                         x_label="time (s)",
                         y_label="angle (degree)",
                         title = "attitude",
                         symbol="r",
                         label="heading",
                         nro_fig=nro_fig + 4)
                    sbes_ek80.plot_xy_add\
                        (ax4,
                         ek80.t_att,
                         ek80.att_att[:,1] * 180. / np.pi,
                         symbol="g",
                         label="pitch")
                    sbes_ek80.plot_xy_add\
                        (ax4,
                         ek80.t_att,
                         ek80.att_att[:,2] * 180. / np.pi,
                         symbol="b",
                         label='roll')
                    ax4.legend()
                    
                        
                if ek80.t_heave is not None:
                    ff5, ax5 = sbes_ek80.plot_xy\
                        (ek80.t_heave,
                         ek80.heave_heave,
                         x_label="time (s)",
                         y_label="heave (m)",
                         title = "heave",
                         symbol="b",
                         nro_fig=nro_fig + 5)
                    
                if ek80.t_pos is not None:
                    
                    ff6, ax6 = sbes_ek80.plot_xy\
                        (ek80.pos_pos[:,0],
                         ek80.pos_pos[:,1],
                         x_label="longitude (°)",
                         y_label="latitude (°)",
                         title = "trajectory",
                         symbol="b",
                         nro_fig=nro_fig + 6)
                    
                    ff7, ax7 = sbes_ek80.plot_xy\
                        (ek80.t_pos,
                         ek80.pos_pos[:,2],
                         x_label="time (s)",
                         y_label="altitude (m)",
                         title = "altitude",
                         symbol="b",
                         nro_fig=nro_fig + 7)

                if ek80.t_course is not None:
                    ff8, ax8 = sbes_ek80.plot_xy\
                        (ek80.t_course,
                         ek80.speed_course,
                         x_label="time (s)",
                         y_label="speed (m)",
                         title = "speed",
                         symbol="b",
                         nro_fig=nro_fig + 8)
                
                    ff9, ax9 = sbes_ek80.plot_xy\
                        (ek80.t_course,
                         ek80.heading_course * 180. / np.pi,
                         x_label="time (s)",
                         y_label="course (m)",
                         title = "course (degrees)",
                         symbol="b",
                         nro_fig=nro_fig + 9)
                

                plt.show()
                #plt.pause(0.1)

            print("batch {:d}: found {:d} pings".\
                  format(nro_batch, ek80.t_ping.size))
                    
            if q_export is not None:
            
                if q_normalize == "TS":
                    label = "TS"
                elif q_normalize == "SV":
                    label = "SV"
                else:
                    label = "raw"

                # Date début/fin de l'enregistrement
                t1 = ek80.reference_time + datetime.timedelta\
                    (seconds=ek80.t_ping[0])
                t2 = ek80.reference_time + datetime.timedelta\
                    (seconds=ek80.t_ping[-1])
                
                file_out = "{:s}_{:s}_{:s}_{:s}{:s}".format\
                    (radix_out, label,
                     t1.strftime("%Y%m%d_%H%M%S"),
                     t2.strftime("%H%M%S"),
                     q_export)
                                                  
                                            
                # Interpolation des positions, pilonnements, attitudes
                # aux dates d'émission
                # Interpolation de la position
                if ek80.t_pos is not None and ek80.t_pos.size > 1: 
                    f = scipy.interpolate.interp1d\
                        (ek80.t_pos,
                         ek80.pos_pos,
                         kind="linear",
                         axis=0,
                         bounds_error=False,
                         fill_value=(ek80.pos_pos[0,:],
                                     ek80.pos_pos[-1,:]),
                         assume_sorted=True)
                    position = f(ek80.t_ping)
                else:
                    position = ek80.pos_pos
                
                # Interpolation de l'attitude
                if ek80.t_att is not None and ek80.t_att.size > 1:
                    f = scipy.interpolate.interp1d\
                        (ek80.t_att,
                         ek80.att_att,
                         kind="linear",
                        axis=0,
                         bounds_error=False,
                         fill_value=(ek80.att_att[0,:],
                                     ek80.att_att[-1,:]),
                         assume_sorted=True)
                    attitude= f(ek80.t_ping)
                else:
                    attitude = ek80.att_att

                if ek80.t_heave is not None and ek80.t_heave.size > 1:
                    # Interpolation du pilonnement
                    f = scipy.interpolate.interp1d\
                        (ek80.t_heave,
                         ek80.heave_heave,
                         kind="linear",
                         axis=0,
                         bounds_error=False,
                         fill_value=(ek80.heave_heave[0],
                                     ek80.heave_heave[-1]),
                         assume_sorted=True)
                    heave= f(ek80.t_ping)
                else:
                    heave = ek80.heave_heave

                if ek80.t_course is not None and ek80.t_course.size > 1:
                    # Interpolation de la vitesse et du cap
                    f = scipy.interpolate.interp1d\
                        (ek80.t_course,
                         ek80.speed_course,
                         kind="linear",
                         axis=0,
                         bounds_error=False,
                         fill_value=(ek80.speed_course[0],
                                     ek80.speed_course[-1]),
                         assume_sorted=True)
                    speed= f(ek80.t_ping)

                    f = scipy.interpolate.interp1d\
                        (ek80.t_course,
                         ek80.heading_course,
                         kind="linear",
                         axis=0,
                         bounds_error=False,
                         fill_value=(ek80.heading_course[0],
                                     ek80.heading_course[-1]),
                         assume_sorted=True)
                    course= f(ek80.t_ping)

                else:
                    speed = ek80.speed_course
                    course = ek80.heading_course

                if q_export == ".mat":
                    # Export en format Matlab
                    date_ = np.empty((6,), np.float32)
                    date_[0] = ek80.reference_time.year
                    date_[1] = ek80.reference_time.month
                    date_[2] = ek80.reference_time.day
                    date_[3] = ek80.reference_time.hour
                    date_[4] = ek80.reference_time.minute
                    date_[5] = ek80.reference_time.second\
                        + ek80.reference_time.microsecond / 1.e6
                    
                    data = {}
                    data["data_type"] = label
                    data["data"] = ek80.amplitude_ping.T
                    data["ref_time"] = date_
                    data["time"] = ek80.t_ping
                    data["range_time"] = ek80.time_range
                    data["position"] = position
                    data["attitude"] = attitude
                    data["heave"] = heave
                    data["svp"] = ek80.sound_velocity_profile
                    data["sound_velocity"] = ek80.sound_velocity
                    data["speed"] = speed
                    data["course"] = course
                    
                    if q_angle == True:
                        data["angle_athwartship"] = ek80.angle_ping[...,0].T
                        data["angle_longitudinal"]\
                            = ek80.angle_ping[...,1].T
                        data["beamwidth_athwartship"]= ek80.beamwidth_athwartship
                        data["beamwidth_alongship"]= ek80.beamwidth_alongship

                        \
                    scipy.io.savemat\
                        (os.path.join(directory_out, file_out), data)


                elif q_export == ".npz":
                    # Export en format Python
                    if q_angle == True:
                        np.savez\
                            (os.path.join(directory_out, file_out),
                             data_type = label,
                             data = ek80.amplitude_ping.T,
                             angle_athwartship = ek80.angle_ping[...,0].T,
                             angle_longitudinal = ek80.angle_ping[...,1].T,
                             beamwidth_athwartship = ek80.beamwidth_athwartship,
                             beamwidth_alongship = ek80.beamwidth_alongship,
                             ref_time = ek80.reference_time,
                             time = ek80.t_ping,
                             range_time = ek80.time_range,
                             position = position,
                             attitude = attitude,
                             heave = heave,
                             course = course,
                             speed = speed,
                             svp = ek80.sound_velocity_profile,
                              sound_velocity = ek80.sound_velocity)
                    else:
                        np.savez\
                            (os.path.join(directory_out, file_out),
                             data_type = label,
                             data = ek80.amplitude_ping.T,
                             ref_time = ek80.reference_time,
                             time = ek80.t_ping,
                             range_time = ek80.time_range,
                             position = position,
                             attitude = attitude,
                             heave = heave,
                             course = course,
                             speed = speed,
                             svp = ek80.sound_velocity_profile,
                             sound_velocity = ek80.sound_velocity)

            if q_iterate == False or q_eof == True:
                # On arrête là (on passe au fichier d'entrée suivant)
                break
            else:
                # On prépare la sortie suivante
                ping0_ = ek80.n0_ping
                nro_batch += 1
                    
