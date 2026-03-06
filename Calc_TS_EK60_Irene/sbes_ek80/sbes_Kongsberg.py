import os
import struct
import datetime
import math
import xml.etree.ElementTree as ET

import numpy as np

class EK60:

    def __init__(self, directory, filelist, ref_time=None):

        self.q_angle = None
        
        self.directory = directory
        self.filelist = filelist
        
        self.type_set = set()
        self.nmea_set = set()

        if isinstance(filelist, (tuple, list)) == False:
            filelist = (filelist,)
                     
        # Ouverture des fichiers
        self.f_list = []
        self.time_ = {}
        self.ref_time = ref_time

        for filename in filelist:
            f = open( os.path.join(directory, filename), 'rb' )
            self.f_list.append(f)            

        self.f_number = 0
        self.f = self.f_list[self.f_number]

    def rewind(self):
        for f in self.f_list:
            if f is not None:
                f.seek(0, os.SEEK_SET)
        self.f_number = 0

    def close(self):
        for i_ in range( len( self.f_list ) ):
            if self.f_list[ i_ ] is not None:
                self.f_list[ i_ ].close()
                self.f_list[ i_ ] = None
        self.f_number = 0
        self.f = None

    def windows_time_to_ref_time(self, ref_time=None):
        """ Conversion d'une date windows vers une date et un temps UTC
          windows time correspond à 2 entiers 32 bits : 
          le premier est le LSB et le second le MSB
          C'est le temps écoulé depuis le premier janvier 1601 à 0h a une 
        résolution de 100 ns
        Il retourne la date et le temps UTC sous forme datetime Python """
        # Conversion du temps windows au temps Unix
        # Il y a 134774 jours ou 11 644 473 600 secondes
        unix_time= (self.time_["high date time"] * (2 ** 32)\
                    + self.time_["low date time"]) * 1e-7\
                    -11644473600

        if ref_time is None:
            return datetime.datetime.utcfromtimestamp(unix_time)
        else:
            return (datetime.datetime.utcfromtimestamp(unix_time)
                    - ref_time).total_seconds()
        
    def read_packet(self, wished_packet=None):

        try:
            n_dgm, = struct.unpack('<I',self.f.read(4)) #Longueur du paquet
        except (struct.error, EOFError):
            # Le paquet est fini, on regarde si il y a un suivant
            self.f_number += 1
            if self.f_number == len( self.f_list ):
                raise EOFError
            self.f = self.f_list[ self.f_number ]
            try:
                n_dgm, = struct.unpack('<I',self.f.read(4)) #Longueur du paquet
            except (struct.error, EOFError):
                raise EOFError
        
        # Lecture du datagramme
        type_, self.time_["low date time"], self.time_["high date time"]\
            = struct.unpack('<4s2I',self.f.read(12))
        #print(type_)

        
        n = 12
        if self.ref_time is None:
            self.ref_time = self.windows_time_to_ref_time()
            self.time = 0.
        else:
            self.time = self.windows_time_to_ref_time(self.ref_time)

        if type_ not in self.type_set:
            self.type_set.add(type_)
            print("packet type:", type_)            
        if wished_packet is None or type_ in wished_packet:

            if type_ == b"CON0":
                self.q_angle = []
                self.config = {}
                #self.config["sounder type"]
                self.config["survey name"],\
                    self.config["transect name"],\
                    self.config["sounder name"],\
                    self.config["motion x"],\
                    self.config["motion y"],\
                    self.config["motion z"],\
                    self.config["spare 1"],\
                    self.config["transducer count"],\
                    = struct.unpack("<128s128s128s3f116sI",
                                    self.f.read(516))
                n += 516
                self.config["transducer"] = []
                for i_ in range(self.config["transducer count"]):
                    x = {}
                    x["channel id"],\
                        x["beam type"],\
                        x["frequency"],\
                        x["gain"],\
                        x["equivalent beam angle"],\
                        x["beamwidth along ship"],\
                        x["beamwidth athwart ship"],\
                        x["angle sensitivity along ship"],\
                        x["angle sensitivity athwart ship"],\
                        x["angle offset along ship"],\
                        x["angle offset athwart ship"],\
                        x["pos x"],\
                        x["pos y"],\
                        x["pos z"],\
                        x["dir x"],\
                        x["dir y"],\
                        x["dir z"],\
                        = struct.unpack("<128sI15f", self.f.read(192))
                    x["pulse length table"] = np.fromfile(self.f, "<f", 5)
                    x["spare 2"] = struct.unpack("<8s", self.f.read(8))
                    x["gain table"] = np.fromfile(self.f, "<f", 5)
                    x["spare 3"] = struct.unpack("<8s", self.f.read(8))
                    x["sa correction table"] = np.fromfile(self.f, "<f", 5)
                    x["spare 4"] = struct.unpack("<52s", self.f.read(52))
                    n += 320
                    self.config["transducer"].append(x)
                    if x["beam type"] == 1:
                        self.q_angle.append(True)
                    else:
                        self.q_angle.append(False)

                self.f.seek(n_dgm - n + 4, os.SEEK_CUR)

            elif type_ == b"NME0":

                n_nmea = n_dgm - n
                nmea_ = self.f.read(n_nmea)

                self.nmea = nmea_.decode().rstrip("\x00") 
                self.f.seek(4, os.SEEK_CUR)

                if self.nmea[:6] not in self.nmea_set:
                    self.nmea_set.add(self.nmea[:6])

                return type_

            elif type_ == b"RAW0":
                self.raw = {}
                self.raw["channel"],\
                self.raw["mode"],\
                self.raw["transducer depth"],\
                self.raw["frequency"],\
                self.raw["transmit power"],\
                self.raw["pulse length"],\
                self.raw["bandwidth"],\
                self.raw["sample interval"],\
                self.raw["sound velocity"],\
                self.raw["absorption coefficient"],\
                self.raw["heave"],\
                self.raw["tx roll"],\
                self.raw["tx pitch"],\
                self.raw["temperature"],\
                self.raw["spare 1"],\
                self.raw["spare 2"],\
                self.raw["rx roll"],\
                self.raw["rx pitch"],\
                self.raw["offset"],\
                self.raw["count"],\
                = struct.unpack("<2H12f2H2f2I", self.f.read(72))
                n += 72

                nro_channel = self.raw["channel"] - 1
                ns = self.raw["count"]
                power = np.fromfile(self.f, "<i2", ns)\
                          .astype("f4")
                self.raw["power"] = power * 10. * math.log10(2) / 256
                n += ns * 2
                
                if self.q_angle is not None\
                   and self.q_angle[nro_channel] == True:
                    angle_ = np.fromfile(self.f, "<b", 2 * ns)\
                        .astype("f4")
                    self.raw["phase"] = angle_ * 180 / 128
                    self.raw["phase"].shape = (-1, 2)
                    n += ns * 2
                else:
                    self.raw["phase"] = None
                self.f.seek(n_dgm - n + 4, os.SEEK_CUR)

            elif type_ == b'DEP0':
                # Fond détecté
                # ------------
                self.depth = []
            
                n_channel, = struct.unpack('I',self.f.read(4))
                for i_ in range(n_channel):
                    depth_ = {}
                    depth_["depth"],\
                        depth_["bs"],\
                        depth_["extra"],\
                        = struct.unpack('3f',self.f.read(12))
                    self.depth.append(depth_)
                n += 4 + n_channel * 12
                self.f.seek(n_dgm - n + 4, os.SEEK_CUR) 
                
            elif type_ == b'TAG0':
                # Annotation
                # -----------
                n_size = n_dgm - 8
                self.tag, = struct.unpack('{:d}s'.format(n_size),
                                          self.f.read(n_size))

            elif type_ == b'SVP0':
                # Bathycélérité entrée
                #----------------------
                n_size = n_dgm - 8
                (self.svp,)=struct.unpack('{:d}s'.format(n_size),
                                          self.f.read(n_size))

            else:
                print("type", type_)
                print("unknown type {:s}".format(type_.decode()))
                self.f.seek(n_dgm - n + 4, os.SEEK_CUR)
            return type_
                      
        else:
            self.f.seek(n_dgm - n + 4, os.SEEK_CUR)

class EK80:

    def __init__(self, directory, filelist, ref_time=None, q_debug=None):
        """
        q_debug = None : on traite par défaut : lecture sous forme d'arbre des fichiers xml et 
              extraction de quelques infos
              = "raw" : affichage des fichiers bruts xml pendant la lecture
              = "tree" : affichage sous forme de dictionnaires des données xml
        """

        self.q_xml_raw = False
        self.q_xml_print = False
        self.q_xml_set = False
        self.q_nmea = False
        self.q_nmea_set = False

        if q_debug is not None:
            if "xml raw" in q_debug:
                self.q_xml_raw = True
            if "xml" in q_debug:
                self.q_xml_print = True
            if "xml set" in q_debug:
                self.q_xml_set = True
            if "nmea" in q_debug:
                self.q_nmea = True
            if "nmea set" in q_debug:
                self.q_nmea_set = True
            
        self.directory = directory
        self.filelist = filelist
        
        self.type_set = set()
        self.nmea_set = set()
        self.xml_set = set()

        if isinstance(filelist, (tuple, list)) == False:
            filelist = (filelist,)
                     
        # Ouverture des fichiers
        self.f_list = []
        self.time_ = {}
        self.ref_time = ref_time

        for filename in filelist:
            f = open( os.path.join( directory, filename), 'rb' )
            self.f_list.append(f)            

        self.f_number = 0
        self.f = self.f_list[self.f_number]

    def rewind(self):
        for f in self.f_list:
            if f is not None:
                f.seek(0, os.SEEK_SET)
        self.f_number = 0

    def close(self):
        for i_ in range( len( self.f_list ) ):
            if self.f_list[ i_ ] is not None:
                self.f_list[ i_ ].close()
                self.f_list[ i_ ] = None
        self.f_number = 0
        self.f = None

    def windows_time_to_ref_time(self, ref_time=None):
        """ Conversion d'une date windows vers une date et un temps UTC
          windows time correspond à 2 entiers 32 bits : 
          le premier est le LSB et le second le MSB
          C'est le temps écoulé depuis le premier janvier 1601 à 0h a une 
        résolution de 100 ns
        Il retourne la date et le temps UTC sous forme datetime Python """
        # Conversion du temps windows au temps Unix
        # Il y a 134774 jours ou 11 644 473 600 secondes
        unix_time= (self.time_["high date time"] * (2 ** 32)\
                    + self.time_["low date time"]) * 1e-7\
                    -11644473600

        if ref_time is None:
            return datetime.datetime.utcfromtimestamp(unix_time)
        else:
            return (datetime.datetime.utcfromtimestamp(unix_time)
                    - ref_time).total_seconds()
        
    def read_packet(self, wished_packet=None):

        try:
            n_dgm, = struct.unpack('<I',self.f.read(4)) #Longueur du paquet
        except (struct.error, EOFError):
            # Le paquet est fini, on regarde si il y a un suivant
            self.f_number += 1
            if self.f_number == len( self.f_list ):
                raise EOFError
            self.f = self.f_list[ self.f_number ]
            try:
                n_dgm, = struct.unpack('<I',self.f.read(4)) #Longueur du paquet
            except (struct.error, EOFError):
                raise EOFError
        
        # Lecture du datagramme
        type_, self.time_["low date time"], self.time_["high date time"]\
            = struct.unpack('<4s2I',self.f.read(12))
        #print(type_)

        
        n = 12
        if self.ref_time is None:
            self.ref_time = self.windows_time_to_ref_time()
            self.time = 0.
        else:
            self.time = self.windows_time_to_ref_time(self.ref_time)

        if type_ not in self.type_set:
            self.type_set.add(type_)
            #print("packet type:", type_)            
        if wished_packet is None or type_ in wished_packet:

            if  type_ == b"TAG0":
                # Annotation texte
                self.txt = self.f.read(n_dgm - 12).decode()
                n = n_dgm

                self.f.seek(4, os.SEEK_CUR)

            elif  type_ == b"XML0":
                #print(n_dgm)
                self.xml_raw = self.f.read(n_dgm-12).decode()
                n = n_dgm

                if self.q_xml_raw == True:
                    # Affichage des données brutes
                    print("#-------------------------------------------------")
                    print(self.xml_raw)
                    print("#-------------------------------------------------")

                # Transformation en arbre des données xml
                xml = ET.fromstring(self.xml_raw.rstrip("\x00"))
                self.xml_type = xml.tag

                if self.q_xml_print == True:
                    print("#----------------------------{:s}".format(xml.tag))

                self.tab = ""
                self.xml_tree = {}
                self.xml_to_tree(self.xml_tree, xml)

                if self.xml_type == "Configuration":
                    self.configuration = self.xml_tree

                    if "configuration" not in self.xml_set:
                        self.xml_set.add("configuration")
                        if self.q_xml_set == True:
                            self.xml_to_print(xml)
                    
                elif self.xml_type == "InitialParameter":
                    self.initial_parameter = self.xml_tree

                    if "initial parameter" not in self.xml_set:
                        self.xml_set.add("initial parameter")
                        if self.q_xml_set == True:
                            self.xml_to_print(xml)

                    
                elif self.xml_type == "Environment":
                    self.environment = self.xml_tree

                    if "environment" not in self.xml_set:
                        self.xml_set.add("environment")
                        if self.q_xml_set == True:
                            self.xml_to_print(xml)

                    
                elif self.xml_type == "Sensor":
                    self.sensor = self.xml_tree

                    if "sensor" not in self.xml_set:
                        self.xml_set.add("sensor")
                        if self.q_xml_set == True:
                            self.xml_to_print(xml)
                    
                elif self.xml_type == "PingSequence":
                    self.ping_sequence = self.xml_tree

                    if "ping sequence" not in self.xml_set:
                        self.xml_set.add("ping sequence")
                        if self.q_xml_set == True:
                            self.xml_to_print(xml)

                    
                elif self.xml_type == "Parameter":
                    self.parameter = self.xml_tree
                    
                    if "parameter" not in self.xml_set:
                        self.xml_set.add("parameter")
                        if self.q_xml_set == True:
                            self.xml_to_print(xml)

                    
                else:
                    print("Unknow XML type: {:s}".format(self.xml_type))
                    
                self.f.seek(4, os.SEEK_CUR)

            elif  type_ == b"MRU0":
                self.mru = {}
                self.mru["heave"],\
                self.mru["roll"],\
                self.mru["pitch"],\
                self.mru["heading"],\
                = struct.unpack("4f", self.f.read(16))
                self.mru["heading"] *= np.pi / 180.
                self.mru["pitch"] *= np.pi / 180.
                self.mru["roll"] *= np.pi / 180.
                
                n += 16
                                
                self.f.seek(n_dgm - n + 4, os.SEEK_CUR)

            elif type_ == b"MRU1":
                print("MRU1 not yet written")
                self.f.seek(n_dgm - n + 4, os.SEEK_CUR)
                
            elif type_ == b"NME0":

                n_nmea = n_dgm - n
                nmea_ = self.f.read(n_nmea)
                self.nmea = nmea_.decode().rstrip("\x00")
                
                if self.q_nmea == True:
                    print(self.nmea)
                
                self.f.seek(4, os.SEEK_CUR)

                if self.nmea[:6] not in self.nmea_set:
                    self.nmea_set.add(self.nmea[:6])
                    if self.q_nmea_set == True:
                        print("nmea type:", self.nmea[:6])

            elif type_ == b"RAW3":
 
                self.raw = {}
                self.raw["channel id"],\
                    self.raw["data type"],\
                    self.raw["offset"],\
                    ns,\
                    = struct.unpack("128sh2x2I", self.f.read(140))
                n += 140

                q_power = self.raw["data type"] & 1
                q_angle = (self.raw["data type"] >> 1) & 1
                q_cmplx_16 = (self.raw["data type"] >> 2) & 1
                q_cmplx_32 = (self.raw["data type"] >> 3) & 1
                n_values = int(self.raw["data type"] >> 8)

                if q_power == 1:
                    # Conversion documentée par Kongsberg
                    # ek80_interface_en_a4.pdf
                    power = np.fromfile(self.f, "<i2", ns).astype("f4")
                    self.raw["power"] = power * 10. * math.log10(2) / 256 
                    n += ns * 2
                    
                if q_angle == 1:
                    # Conversion documentée par Kongsberg
                    # ek80_interface_en_a4.pdf
                    # Reste une partie à compléter avec le paramètre
                    # angle sensitivy 
                    angle_ = np.fromfile(self.f, "<b", 2 * ns).astype("f4")
                    self.raw["phase"] = angle_ * 180 / 128
                    self.raw["phase"].shape = (-1, 2)
                    n += ns *2
                    
                if q_cmplx_32 == 1:
                    
                    self.raw["data"] = np.fromfile\
                        (self.f, np.complex64, n_values * ns)
                    self.raw["data"].shape = (self.raw["count"], -1)
                    n += 8 * n_values * ns

                elif q_cmplx_16 == 1:
                    x = np.fromfile\
                        (self.f, np.float16, n_values * ns * 2)
                    x.shape = (ns, n_values, 2)
                    self.raw["data"] = x[...,0] + 1j * x[...,1]
                    n += 4 * n_values * ns
                    
                                        
                    n += 4 * 2 * self.raw["count"]
                    
                #assert self.dgm == 4
                self.f.seek(n_dgm -n + 4, os.SEEK_CUR)

                
            elif type_ == b'TAG0':
                # Annotation
                # -----------
                n_size = n_dgm - 8
                self.tag, = struct.unpack('{:d}s'.format(n_size),
                                          self.f.read(n_size))

            elif type_ == b'SVP0':
                # Bathycélérité entrée
                #----------------------
                n_size = n_dgm - 8
                (self.svp,)=struct.unpack('{:d}s'.format(n_size),
                                          self.f.read(n_size))

            else:
                print("type", type_)
                print("unknown type {:s}".format(type_.decode()))
                self.f.seek(n_dgm - n + 4, os.SEEK_CUR)
            return type_
                      
        else:
            self.f.seek(n_dgm - n + 4, os.SEEK_CUR)

    def xml_to_tree(self, dict_, xml):
        if self.q_xml_print==True:
            print("{:s}{:s}:".format(self.tab,xml.tag))
            self.tab += "\t"
            
        for k, v in xml.attrib.items():
            if self.q_xml_print == True:
                print("{:s}{:s}: {:s}".format(self.tab, k, v))
            dict_[k] = v

        #dict_["node"] = []
        for i_, xml_ in enumerate(xml):
            #if self.q_xml_print == True:
            #    print("{:s}{:s}:".format(self.tab,xml_.tag))
            #self.tab += "\t"        

            #u = {}
            #self.xml_to_tree(u, xml_)
            #dict_["node"].append(u)
            if xml_.tag not in dict_.keys():
                dict_[xml_.tag] = {}
                self.xml_to_tree(dict_[xml_.tag], xml_)
            else:
                if type(dict_[xml_.tag]) != tuple:
                    dict_[xml_.tag] = [ dict_[xml_.tag] ]
                dict_[xml_.tag].append({})
                self.xml_to_tree(dict_[xml_.tag][-1], xml_)
        if self.q_xml_print==True:
            self.tab = self.tab[:-1]

    def xml_to_print(self, xml):
        print("{:s}{:s}:".format(self.tab,xml.tag))
        self.tab += "\t"
        
        for k, v in xml.attrib.items():
            print("{:s}{:s}: {:s}".format(self.tab, k, v))
            
        for i_, xml_ in enumerate(xml):
            #print("{:s}{:s}:".format(self.tab,xml_.tag))
            #self.tab += "\t"        
            self.xml_to_print(xml_)
        self.tab = self.tab[:-1]

            
if 0:        
    def decode_xml(self):
        """ Décode les différents type de xml """
        # Décodage du xml
        xml = ET.fromstring(self.xml.rstrip("\x00"))
        xml_type = xml.tag

        if xml.tag == "Configuration":
            
            self.configuration = {}
            # Configuration Attribs
            for k, v in xml.attrib.items():
                print("extra attributes in Configuration: {:s}".format(self.k))

            # Configuration childs
            for i_,x0 in enumerate(xml):
                
                if x0.tag == "Header":
                    # Configuration/header Attribs
                    for k, v in x0.attrib.items():
                        if k == "Copyright":
                            self.configuration["copyright"] = v
                        elif k == "ApplicationName":
                            self.configuration["application name"] = v
                        elif k == "Version":
                            self.configuration["version"] = v
                        elif k == "FileFormatVersion":
                            self.configuration["file format version"] = v
                        elif k == "TimeBias":
                            self.configuration["time bias"] = float(v)
                        else:
                            print("extra attributes in Configuration/Header: {:s}".format(k))
                            
                    # Configuration/header childs                            
                    for i_,x1 in enumerate(x0):
                        print('extra childs in configuration"\
                        +"/Header{:s}"'.format(x1.tag))
                                
                elif x0.tag == "ActivePingMode":
                    # Configuration/activePingMode attribs
                    for k, v in x0.attrib.items():
                        if k == "Mode":
                            self.configuration["mode"] = v
                        else:
                            print("extra attributes in Configuration/ActivePingMode: {:s}".format(k))

                    # Configuration/activePingMode childs
                    for i_,x1 in enumerate(x0):
                        print('extra childs in configuration/ActivePingMode:{:s}"'.format(x1.tag))

                elif x0.tag == "Transceivers":
                    # Configuration/Transceivers attribs
                    for k, v in x0.attrib.items():
                        if k == "MergeOperation":
                            pass
                        else:
                            print("extra attributes in Configuration/Transceivers: {:s}".format(k))

                            
                    # Configuration/Transceivers child
                    self.configuration["transceivers"] = []
                    for i_,x1 in enumerate(x0):
                        
                        if x1.tag == "Transceiver":
                            # Configuration/Transceivers/transceiver attrib
                            u1 = {}
                            for k, v in x1.attrib.items():
                                if k == "TransceiverName":
                                    u1["transceiver name"] = v
                                elif k == "IPAddress":
                                    u1["ip address"] = v
                                elif k == "MarketSegment":
                                    u1["market segment"] = v
                                elif k == "SerialNumber":
                                    u1["serial no"] = v
                                elif k == "Impedance":
                                    u1["impedance"] = float(v)
                                elif k == "Multiplexing":
                                    u1["multiplexing"] = v
                                elif k == "RxSampleFrequency":
                                    u1["rx sample frequency"] = float(v)
                                elif k == "EthernetAddress":
                                    u1["ethernet address"] = v
                                elif k == "Version":
                                    u1["version"] = v
                                elif k == "TransceiverSoftwareVersion":
                                    u1["transceiver software version"] = v
                                elif k == "TransceiverNumber":
                                    u1["transceiver number"] = int(v)
                                elif k == "TransceiverType":
                                    u1["transceiver type"] = v
                                else:
                                    print("extra attributes in Configuration/Transceivers/transceiver: {:s}".format(k))
                                    
                                # Configuration/Transceivers/transceiver child
                                for i_,x2 in enumerate(x1):
                                    u1["channels"] = []

                                    # Configuration/Transceivers/transceiver/channels attrib
                                    for k, v in x2.attrib.items():
                                        print("extra attributes in Configuration/Transceivers/transceiver/channels: {:s}".format(k))

                                    # Configuration/Transceivers/transceiver/channels childs                                        
                                    for i_,x3 in enumerate(x2):

                                        if x3.tag == "channel":
                                            u2 = []
                                            # Configuration/Transceivers/transceiver/channels/channel attrib
                                            for k, v in x3.attrib.items():
                                                if k == "ChannelID":
                                                    u2["channel id"] = v
                                                elif k == "LogicalChannelID":
                                                    u2["logical channel id"] = v
                                                elif k == "ChannelIdShort":
                                                    u2["channel id short"] = v
                                                elif k == "Serial No":
                                                    u2["serial no"] = int(v)
                                                elif k == "MaxTxPowerTransceiver":
                                                    u2["max tx power transceiver"] = float(v)
                                                elif k == "HWChannelConfiguration":
                                                    u2["HW channel configuration"] = float(v)
                                                elif k == "PulseDuration":
                                                    v_ = v.split(";")
                                                    vv = np.empty((len(v_),), np.float32)
                                                    for i_, vv_ in enumerate(v_):
                                                        vv[i_] = float(vv_)
                                                    u2["pulse duration"] = vv
                                                else:
                                                    print("extra attributes in Configuration/Transceivers/transceiver/channels/channel: {:s}".format(k))
                                            # Configuration/Transceivers/transceiver/channels/channel child
                                            for i_,x4 in enumerate(x3):
                                                # Configuration/Transceivers/transceiver/channels/channel/transducer
                                                if x4.tag == "Transducer":
                                                    u3 = {}
                                                    for k, v in x4.attrib.items():
                                                        if k == "TransducerName":
                                                            u3["transducer name"] = v
                                                        elif k == "SerialNumber":
                                                            u3["serial n0"] = int(v)
                                                        elif k == "Frequency":
                                                            u3["frequency"] = float(v)
                                                        elif k == "FrequencyMinimum":
                                                            u3["frequency minimum"] = float(v)
                                                        elif k == "FrequencyMaximum":
                                                            u3["frequency maximum"] = float(v)
                                                        elif k == "MaxTxPowerTransducer":
                                                            u3["max tx power transducer"] = float(v)
                                                        elif k == "Gain":
                                                            v_ = v.split(";")
                                                            vv = np.empty((len(v_),), np.float32)
                                                            for i_, vv_ in enumerate(v_):
                                                                vv[i_] = float(vv_)
                                                            u3["gain"] = vv
                                                        elif k == "SaCorrection":
                                                            v_ = v.split(";")
                                                            vv = np.empty((len(v_),), np.float32)
                                                            for i_, vv_ in enumerate(v_):
                                                                vv[i_] = float(vv_)
                                                            u3["Sa correction"] = vv
                                                        elif k == "EquivalentBeamAngle":
                                                            u3["equivalent beam angle"] = float(v)
                                                        elif k == "DirectivityDropAt2XBeamWidth":
                                                            u3["directivity drop at 2x beamwidth"] = float(v)
                                                        elif k == "AngleSensitivityAlongship":
                                                            u3["angle sensitivity alongship"] = float(v)
                                                        elif k == "AngleSensitivityAthwartship":
                                                            u3["angle sensitivity athwartship"] = float(v)
                                                        elif k == "AngleOffsetAlongship":
                                                            u3["angle offset alongship"] = float(v)
                                                        elif k == "AngleOffsetAthwartship":
                                                            u3["angle offset athwartship"] = float(v)
                                                        elif k == "BeamWidthAlongship":
                                                            u3["beamwidth alongship"] = float(v)
                                                        elif k == "BeamWidthAthwartship":
                                                            u3["beam width athwartship"] = float(v)
                                                        else:
                                                            print("extra attributes in Configuration/Transceivers/transceiver/channels/channel/transceiver: {:s}".format(k))
                                                    # Configuration/Transceivers/transceiver/channels/channel/transducer
                                                    for k, v in x4:                
                                                        print('extra childs in configuration/Transceivers/Transceiver/Channels/Channel/transducer: {:s}"'.format(x4.child))
                                                u2["transducer"] = u3
                                            u1["channels"].append(u2)
                        self.configuration["transceivers"].append(u1)
                                                         
                elif x0.tag == "Transducers":
                    # Configuration/transducers attrib
                    self.configuration["transducers"] = []
                    for k, v in x0.attrib.items():
                        if k == "MergeOperation":
                            pass
                        else:
                            print("extra attributes in Configuration/Transducers: {:s}".format(k))

                    # Configuration/transducers child
                    for i_, x1 in enumerate(x0):
                        if x1.tag == "Transducer":
                            # Configuration/transducers/transducer attrib
                            u1 = {}
                            for k, v in x1.attrib.items():
                                if k == "TransducerName":
                                    u1["transducer name"] = v
                                elif k == "TransducerSerialNumber":
                                    u1["transducer serial number"] = int(v)
                                elif k == "TransducerCustomName":
                                    u1["transducer custom name"] = v
                                elif k == "TransducerMounting": 
                                    u1["transducer mounting"] = v
                                elif k == "TransducerOffsetX":
                                    u1["transducer offset X"] = float(v)
                                elif k == "TransducerOffsetY":
                                    u1["transducer offset Y"] = float(v)
                                elif k == "TransducerOffsetZ":
                                    u1["transducer offset Z"] = float(v)
                                elif k == "TransducerAlphaX":
                                    u1["transducer alpha X"] = float(v)
                                elif k == "TransducerAlphaY":
                                    u1["transducer alpha Y"] = float(v)
                                elif k == "TransducerAlphaZ":
                                    u1["transducer alpha Z"] = float(v)
                                else:
                                    print("extra attributes in Configuration/Transducer/Transducer: {:s}".format(k))
                                    
                            # Configuration/Transducers/Transducer child
                            for i_, x2 in enumerate(x1):                
                                print("extra childs in configuration/Transducers/Transducer: {:s}".format(x2.tag))
                                    
                        self.configuration["transducers"].append(u1)
                                    
                elif x0.tag == "ConfiguredSensors":
                    # Configuration/ConfiguredSensors
                    self.configuration["configured sensors"] = {}
                    for k, v in x0.attrib.items():
                        print("extra attributes in Configuration/Transducers: {:s}".format(k))

                    
                else:
                    print('extra childs in Configuation: {:s}"'.format(x0.tag))                
            
        elif xml.tag == "InitialParameter":
            
            
            pass
        elif xml.tag == "Environment":
            pass
        elif xml.tag == "PingSequence":
            pass
        elif xml.tag == "Sensor":
            pass
        elif xml.tag == "Parameter":
            pass
        else:
            print("Unknown tag parameter: {:s}".format(xml_type))
            
        print("*********")
        # print("xml", xml_type)
        import pdb;pdb.set_trace()        
        
        return xml_type
    
