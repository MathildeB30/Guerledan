#-*- coding: utf-8 -*-
""" Décodage de trames NMEA 0104

WARNING: Message MUST be a STRING

Les trames actuellement interprétées sont :

GPS standards
-------------
xxGGA : position
xxGGK : position
xxGLL : position
xxVTG : vitesse
GPZDA : date
GPGST et GNGST  (Trimble ?) GPS statistiques

Centrales de cap standards
---------------------------
HEHDT : cap

Capteurs FURUNO
-----------------
PFEC (GPatt / GPhve ) 

Centrales Octans
------------------
1er paquet HEHDT
PHTRO : roulis/tangage Octans
PHLIN : surge/sway/heave Octans
PHSPD : vitesses surge/sway/heave Octans
PHCMP : valeurs de compensation
PHINF : statut Octans
PHZDA : time Octans
Les trames Octans sont directement mises dans le repère "hydro"

Centrales Phins
------------------

1er paquet HEHDT
PIXSE, ATITUD :
PIXSE, POSITI :
PIXSE, SPEED_ :
PIXSE, UTMWGS :
PIXSE, HEAVE_ :
PIXSE, STDHRP :
PIXSE, STDPOS :
PIXSE, STDSPD :
PIXSE, TIME__ :
PIXSE, GPSIN_ :
PIXSE, UTCIN_ :
PIXSE, ALGSTS :
PIXSE, STATUS :
PIXSE, HT_STS :

Attitude RDI et PosMV
---------------------
PRDID : attitude

Centrales PosMV
-----------------
PASHR : attitude (Tait Bryan ou TSS)

GPS Trimble ? (Panopée)
-------------------
IIXTE :
IIRMB : 
IIRMC : 
PCMPN : 
IIAPB :



WARNING : on ne gère pas le type fournisseur pour l'instant : toujours GP !!!
"""

import numpy as np
import math

def NMEA_decode(message):
    """ Décodage d'une trame NMEA 0104
    on suppose qu'il n'y a qu'une trame NMEA 0104 
    qui commence avec le premier caractère comme $"""
    frame={}
    # Début de trame
    if message[0] != '$':
        print("le message n'est pas un code NMEA : {:}".format( message ))
        frame["type"]="000"
        return frame

    # Découpage de la suite de la trame (sauf le checksum)
    l=message[1:].split(",")
    # On enlève le checksum le cas échant
    l[-1]=l[-1].split("*")[0]
    # Type de message
    try:
        if l[0][2:]=="GGA":
            # Position et temps
            frame["type"]=l[0][2:]
            frame["origin"] = l[0][:2]
            # UTC time from midnight
            frame["time"]=float(l[1][0:2])*3600.+\
                float(l[1][2:4])*60.+float(l[1][4:])
            # Latitude en radians
            x=float(l[2][0:2])+float(l[2][2:])/60.
            if l[3]=='N':
                frame["latitude"]=x*math.pi/180.
            else:
                frame["latitude"]=-x*math.pi/180.
            # Longitude en radians
            x=float(l[4][0:3])+float(l[4][3:])/60.
            if l[5]=='E':
                frame["longitude"]=x*math.pi/180.
            else:
                frame["longitude"]=-x*math.pi/180.
            # Gps quality indicator
            frame["quality"]=int(l[6])
            # Satellite number
            frame["satellite nbr"]=int(l[7])
            # HDOP
            frame["hdop"]=float(l[8])
            # Altitude
            frame["altitude"]=float(l[9])
            if l[10]!='M':
                print("Pb on altitude unit on GGA {:s}".format(l[10]))
                      
            # Géoide séparation
            if len(l[11]):
                frame["geoid"]=float(l[11])
            if l[12]!='M':
                print("Pb on geoide separation unit on GGA {:s}"\
                              .format(l[12]))
            # La suite n'est pas pertinente dans le cas d'un GPS absolu
            if frame["quality"]>2:
                # Age of differential date
                frame["age"]=float(l[13])
                # Number of ref. station
                frame["station"]=float(l[14][0:4])

        elif l[0][2:]=="GGK":
            # Position et temps mode RTK (Trimble)
            # Position et temps
            frame["type"]=l[0][2:]
            frame["origin"] = l[0][:2]
            # UTC time from midnight
            frame["time"]=float(l[1][0:2])*3600.+\
                float(l[1][2:4])*60.+float(l[1][4:])
            # UTC date
            frame["date"]=int(l[2])            
            # Latitude en radians
            x=float(l[3][0:2])+float(l[3][2:])/60.
            if l[4]=='N':
                frame["latitude"]=x*math.pi/180.
            else:
                frame["latitude"]=-x*math.pi/180.
            # Longitude en radians
            x=float(l[5][0:3])+float(l[5][3:])/60.
            if l[6]=='E':
                frame["longitude"]=x*math.pi/180.
            else:
                frame["longitude"]=-x*math.pi/180.
            # Gps quality indicator
            frame["quality"]=int(l[7])
            # Satellite number
            if len( l[8] ):
                frame["satelliteNbr"]=int(l[8])
            else:
                frame["satelliteNbr"]=-1
            # DOP
            frame["dop"]=float(l[9])
            # Altitude ellipsoid
            assert l[10][:3]=="EHT"
            frame["altitude"]=float(l[10][3:])
            if l[11]!='M':
                print("Pb on altitude unit on GGK {:s}".format(l[10]))
            
        elif l[0][2:]=="GLL":
            # Position et temps
            frame["origin"] = l[0][:2]
            frame["type"]=l[0][2:]

            # Latitude en radians
            x=float(l[1][0:2])+float(l[1][2:])/60.
            if l[2]=='N':
                frame["latitude"]=x*math.pi/180.
            else:
                frame["latitude"]=-x*math.pi/180.
            # Longitude en radians
            x=float(l[3][0:3])+float(l[3][3:])/60.
            if l[4]=='E':
                frame["longitude"]=x*math.pi/180.
            else:
                frame["longitude"]=-x*math.pi/180.
            # Eventuellement le temps du fix
            if len(l)>5:
                # On a aussi le temps
                frame["q time"] = True
                frame["time"] = float(l[5][0:2])*3600.+\
                float(l[5][2:4])*60.+float(l[5][4:])
            else:
                frame["q time"] = False
                
        elif l[0][2:]=="VTG":
            # Lecture d'une trame cap/vitesse
            frame["type"]=l[0][2:]
            frame["origin"]=l[0][:2]
            # Vrai cap
            frame["true heading"]=float(l[1])
            if l[2]!='T':
                print("Pb on definition of true heading {:s}"
                              .format(l[2]))
            # Cap magnétique : peut ne pas être renseigné !!!
            if len(l[3]):
                frame["mag heading"]=float(l[3])
            if l[4]!='M':
                print("Pb on definition of magnetic heading {:s}"
                              .format(l[4]))
            # speed in knots
            frame["knot speed"]=float(l[5])
            if l[6]!='N':
                print("Pb of knot unit {:s}".format(l[6]))
            # speed in meters
            frame["speed"]=float(l[7])*1000./3600.
            if l[8]!='K':
                print("Pb on speed unit {:s}".format(l[8]))
            # mode indicator
            frame["mode"]=l[9][0]

        elif l[0]=="GPZDA":
            # lecture de la date et de la zone horaire
            frame["type"]=l[0]
            # UTC time from midnight
            frame["time"]=float(l[1][0:2])*3600.+\
                float(l[1][2:4])*60.+float(l[1][4:])
            # Date
            frame["date"]=l[4]+l[3]+l[2]
            # Zone (en heures décimales)
            try:  # Pas toujours mentionné
                frame["zone"]=float(l[5])+float(l[6][0:2])/60.
            except:
                frame["zone"]=None

        elif l[0][2:]=="GST":
            # Lecture d'une trame de statistique pseudorange récepteur GPS
            frame["type"]=l[0][2:]
            frame["origin"] = l[0][:2]
            # UTC time from midnight (en secondes)
            frame["time"]=float(l[1][0:2])*3600.+\
                float(l[1][2:4])*60.+float(l[1][4:])
            # Rms pseudo ranges
            frame["sigmaR"]=float(l[2])
            # Ellipsoid axis major RMS
            frame["sigmaA"]=float(l[3])
            # Ellipsoid axis minor
            frame["sigmaB"]=float(l[4])
            # Ellipsoid horizontal angle
            frame["ellipsoidTheta"]=float(l[5])
            # Latitude RMS
            frame["latRMS"]=float(l[6])
            # Longitude RMS
            frame["lonRMS"]=float(l[7])
            # Altitude RMS
            frame["altitudeRMS"]=float(l[8])

        elif l[0][2:]=="XTE":
            # Erreur par rapport à la piste
            frame["type"] = l[0][2:]
            frame["origin"] = l[0][:2]
            # Général warning
            frame["validity"] = l[1]
            # Loran C cycle lock
            frame["Loran C cycle lock"] = l[2]
            # Distance à la piste
            if l[4] == 'L':
                frame["cross track distance"] = float(l[3])
            elif l[4] == 'R':
                frame["cross track distance"] = -float(l[3])
            # Unités
            frame["unit"] = l[5]
            if l[5] == "N":
                frame["cross track distance"] *= 1852.

        elif l[0][2:] == "RMB":
            # Infos de navigation par rapport à un waypoint
            frame["type"] = l[0][2:]
            frame["origin"] = l[0][:2]
            frame["validity"] = l[1]
            if l[3] == 'L':
                frame["cross track distance"] = float(l[2])
            elif l[3] == 'R':
                frame["cross track distance"] = -float(l[2])
            # Numéro du waypoint
            frame["origin waypoint id"] = l[4]
            frame["destination waypoint id"] = l[5]
            # Position du point de destination
            # Latitude en radians
            x=float(l[6][0:2])+float(l[6][2:])/60.
            if l[7]=='N':
                frame["destination latitude"]=x*math.pi/180.
            else:
                frame["destination latitude"]=-x*math.pi/180.
            # Longitude en radians
            x=float(l[8][0:3])+float(l[8][3:])/60.
            if l[9]=='E':
                frame["destination longitude"]=x*math.pi/180.
            else:
                frame["destination longitude"]=-x*math.pi/180.
            # Conversion en mètre de la distance
            frame["destination range"] = float(l[10])*1852.
            frame["destination bearing"] = float(l[11])*math.pi/180.
            # Vitesse radiale
            frame["closing velocity"] = float(l[12])*1850./3600.
            frame["arrival status"] = l[13] 

        elif l[0][2:] == "RMC":
            # Informations de transit
            frame["type"] = l[0][2:]
            frame["origin"] = l[0][:2]

            # UTC time from midnight
            frame["time"]=float(l[1][0:2])*3600.+\
                float(l[1][2:4])*60.+float(l[1][4:])
            frame["validity"] = l[2]
            # Latitude / longitude
            x=float(l[3][0:2])+float(l[3][2:])/60.
            if l[4]=='N':
                frame["latitude"]=x*math.pi/180.
            else:
                frame["latitude"]=-x*math.pi/180.
            # Longitude en radians
            x=float(l[5][0:3])+float(l[5][3:])/60.
            if l[6]=='E':
                frame["longitude"]=x*math.pi/180.
            else:
                frame["longitude"]=-x*math.pi/180.
            frame["speed"] = float(l[7])*1850./3600.
            frame["heading"] = float(l[8])*math.pi/180.
            frame["date"]=l[9]

            if l[11] == 'E':
                frame["magnetic declination"] = float(l[10])*math.pi/180.
            else:
                frame["magnetic declination"] = -float(l[10])*math.pi/180.

        elif l[0][2:] == "APB":
            # Autopilot
            frame["type"] = l[0][2:]
            frame["origin"] = l[0][:2]
            # Général warning
            frame["validity"] = l[1]
            # Loran C cycle lock
            frame["Loran C cycle lock"] = l[2]
            # Distance à la piste
            if l[4] == 'L':
                frame["cross track distance"] = float(l[3])
            elif l[4] == 'R':
                frame["cross track distance"] = -float(l[3])
            # Unités
            frame["unit"] = l[5]
            if l[5] == "N":
                frame["cross track distance"] *= 1852.
            frame["arrival circle entered"] = l[6]
            frame["perpendicular passed at waypoint"] = l[7]
            frame["origin to destination bearing"] = float(l[8])*math.pi/180.
            frame["origin to destination bearing type"] = l[9]
            frame["destination waypoint id"] = l[10]
            frame["present to destination bearing"] = float(l[11])*math.pi/180.
            frame["present to destination bearing type"] = l[12]
            frame["heading to destination heading"] = float(l[13])*math.pi/180.
            frame["heading to destination heading type"] = l[14]
            
            
        elif l[0]=="PFEC":
            # Matériel FURUNO
            frame["type"]=l[0]
            frame["subType"]=l[1]
            if l[1]=="GPatt":
                # Attitude gps Furuno
                frame["heading"]=float(l[2])*math.pi/180.
                frame["roll"]=float(l[3])*math.pi/180.
                frame["pitch"]=float(l[4])*math.pi/180.

            elif l[1]=="GPhve":
                # Pilonnement gps Furuno ($PFEC,GPhve) 
                frame["heave"]=float(l[2])


        elif l[0]=="PHZDA":
            # lecture de la date et de la zone horaire Octans
            frame["type"]=l[0]
            # UTC time from midnight
            frame["time"]=float(l[1][0:2])*3600.+\
                float(l[1][2:4])*60.+float(l[1][4:])
            # Date
            frame["date"]=int(l[4]+l[3]+l[2])
            # Zone (en heures décimales)
            try:  # Pas toujours mentionné
                frame["zone"]=float(l[5])+float(l[6][0:2])/60.
            except:
                frame["zone"]=None

        elif l[0][2:]=="HDT":
            # lecture du cap
            frame["type"]=l[0]
            # Valeur du cap (conversion en radians)
            frame["heading"]=float(l[1])*math.pi/180.
            # Type de cap
            frame["headingType"]=l[2][0]

        elif l[0]=="PHTRO":
            # lecture du roulis/tangage Octans
            # Conventions : pitch+ bow up
            frame["type"]=l[0]
            # Valeur du pitch
            if l[2]=='P':
                frame["pitch"]=-float(l[1])*math.pi/180.
            elif l[2]=='M':
                frame["pitch"]=float(l[1])*math.pi/180.
            # Valeur du roulis
            if l[4]=='B':
                frame["roll"]=-float(l[3])*math.pi/180.
            elif l[4]=='T':
                frame["roll"]=float(l[3])*math.pi/180.

        elif l[0]=="PRDID":
            # lecture du roulis/tangage RDI
            # Conventions : pitch+ bow up
            frame["type"]=l[0]
            # Valeur du pitch
            frame["pitch"] = float(l[1])*math.pi/180.
            # Valeur du roulis
            frame["roll"]=float(l[2])*math.pi/180.
            # Valeur du cap
            frame["heading"]=float(l[3])*math.pi/180.

        elif l[0]=="PHLIN":
            # Surge/sway/heave
            frame["type"]=l[0]
            frame["surge"]=float(l[1])
            frame["sway"]=-float(l[2])   #Octans : défaut positive port
            frame["heave"]=-float(l[3])  #Octans : défaut positive up

        elif l[0]=="PHSPD":
            # Surge/sway/heave
            frame["type"]=l[0]
            frame["vSurge"]=float(l[1])
            frame["vSway"]=-float(l[2])
            frame["vHeave"]=-float(l[3])

        elif l[0]=="PHCMP":
            # Valeurs de compensation entrées
            frame["type"]=l[0]
            # Latitude en radians
            x=float(l[1][0:2])+float(l[1][2:])/60.
            if l[2]=="N":
                frame["latitude"]=math.pi/180.*x
            else:
                frame["latitude"]=-math.pi/180.*x
            frame["speed"]=float(l[3])

        elif l[0]=="PHINF":
            frame["type"]=l[0]
            frame["status"]=l[1]

        elif l[0]=="PASHR":
            # lecture du roulis/tangage Octans
            # Conventions : pitch+ bow up
            frame["type"]=l[0]
            # UTC time from midnight
            frame["time"]=float(l[1][0:2])*3600.+\
                float(l[1][2:4])*60.+float(l[1][4:])
            # Valeur du cap (en radians)
            frame["heading"] = math.pi/180.*float(l[2])

            assert l[3] == "T"

            # Valeur du roulis (convention Tate Bryant)
            frame["roll"] = math.pi/180.*float(l[4])
            
            # Valeur du pitch (convention Tate Bryant)
            frame["pitch"] = math.pi/180.*float(l[5])

            # Valeur du pilonnement
            frame["heave"] = float(l[6])
            
            # Précision du roulis
            frame["roll accuracy"] = float(l[7])
            # Précision du roulis
            frame["pitch accuracy"] = float(l[8])
            # Précision du roulis
            frame["heading accuracy"] = float(l[9])
            # Type de mesure du heading
            frame["heading mode"] = int(l[10])
            # Validité de l'IMU
            frame["IMU status"] = int(l[11])

        elif l[0]=="PIXSE":
            frame["type"] = l[0]
            # Données d'une PHINS
            frame["subtype"] = l[1]
            if l[1] == "ATITUD":
                frame["roll"] = float(l[2])*math.pi/180.
                frame["pitch"] = float(l[3])*math.pi/180.
            elif l[1] == "POSITI":
                frame["latitude"] = float(l[2])*math.pi/180.
                frame["longitude"] = float(l[3])*math.pi/180.
                frame["altitude"] = float(l[4])
            elif l[1] == "SPEED_":
                frame["east speed"] = float(l[2])
                frame["north speed"] = float(l[3])
                frame["up speed"] = float(l[4])
            elif l[1] == "UTMWGS":
                frame["latitude zone"] = l[2]
                frame["longitude zone"] = int(l[3])
                frame["easting"] = float(l[4])
                frame["northing"] = float(l[5])
                frame["altitude"] = float(l[6])
            elif l[1] == "STDHRP":
                frame["std heading"] = float(l[2])*math.pi/180.
                frame["std pitch"] = float(l[3])*math.pi/180.
                frame["std roll"] = float(l[4])*math.pi/180.
            elif l[1] == "STDPOS":
                frame["std latitude"] = float(l[2])
                frame["std longitude"] = float(l[3])
                frame["std altitude"] = float(l[4])
            elif l[1] == "STDSPD":
                frame["std north speed"] = float(l[2])
                frame["std east speed"] = float(l[3])
                frame["std up speed"] = float(l[4])
            elif l[1] == "HEAVE_":
                frame["surge"] = float(l[2])
                frame["sway"] = float(l[3])
                frame["heave"] = float(l[4])
            elif l[1] == "TIME__":
                frame["time"] = float(l[2][:2])*3600. + float(l[2][2:4])*60.\
                                +float(l[2][4:])
            elif l[1] == "UTCIN_":
                frame["time"] = float(l[2][:2])*3600. + float(l[2][2:4])*60.\
                                +float(l[2][4:])
            elif l[1] == "GPSIN_":
                frame["latitude"] = float(l[2])*math.pi/180.
                frame["longitude"] = float(l[3])*math.pi/180.
                frame["altitude"] = float(l[4])
                frame["time"] = float(l[5][:2])*3600. + float(l[5][2:4])*60.\
                                +float(l[5][4:])
                frame["quality"] = int(l[6])
            elif l[1] == "ALGSTS":
                frame["status 1"] = int(l[2],16)
                frame["status 2"] = int(l[3],16)
            elif l[1] == "STATUS":
                frame["status 1"] = int(l[2],16)
                frame["status 2"] = int(l[3],16)
            elif l[1] == "HT_STS":
                frame["status"] = int(l[2],16)
        else:
            #print("unknown NMEA type {:s}".format(l[0]))
            frame["type"]='000'
            
    except:
        print("message: {:s}".format(message))
        print("Bad message form")
        frame["type"]='000'

    return frame
