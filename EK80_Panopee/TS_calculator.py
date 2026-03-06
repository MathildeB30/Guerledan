import numpy as np
import echopype as ep
import xml.etree.ElementTree as ET
from scipy.signal import find_peaks, peak_widths
import pandas as pd
from pyproj import Transformer

class EK80Processor:
    def __init__(self, calib_file, channel):
        self.calib_file = calib_file
        self.channel = channel
        self.f_s = 1.5e6
        
        # Paramètres de calibration fixes (XML)
        self.z_td_e = self._extract_calib_param('./Calibration/Common/PreviousModelParameters', 'Impedance')[0]
        self.frequencies = self._extract_calib_param('./Calibration/CalibrationResults', 'Frequency')
        self.gain = self._extract_calib_param('./Calibration/CalibrationResults', 'Gain')
        self.n_f_points = len(self.frequencies)
        self.angle_offset_alongship_fn = self._extract_calib_param('./Calibration/CalibrationResults','AngleOffsetAlongship')
        self.angle_offset_athwartship_fn = self._extract_calib_param('./Calibration/CalibrationResults','AngleOffsetAthwartship')
        self.beam_width_alongship_fn = self._extract_calib_param('./Calibration/CalibrationResults','BeamWidthAlongship')
        self.beam_width_athwartship_fn = self._extract_calib_param('./Calibration/CalibrationResults','BeamWidthAthwartship')
        
        # Stockage des données du fichier courant
        self.ed = None
        self.groups = []

    def _extract_calib_param(self, xpath, tag):
        tree = ET.parse(self.calib_file)
        root = tree.getroot()
        params_node = root.find(xpath)
        tag_element = params_node.find(tag)
        raw_text = tag_element.text
        if raw_text is None or not raw_text.strip(): return []
        if ';' in raw_text:
            return [float(val) for val in raw_text.split(';') if val.strip()]
        try: return float(raw_text.strip())
        except: return raw_text.strip()

    def _trim_nans(self, data, axis):
        other_axes = tuple(i for i in range(data.ndim) if i != axis)
        is_valid = np.any(~np.isnan(data), axis=other_axes)
        valid_indices = np.where(is_valid)[0]
        if valid_indices.size > 0:
            slc = [slice(None)] * data.ndim
            slc[axis] = slice(0, valid_indices[-1] + 1)
            return data[tuple(slc)]
        return data

    def load_raw(self, raw_path):
        """Charge le fichier et prépare les variables d'environnement."""
        self.ed = ep.open_raw(raw_path, sonar_model="EK80")
        bg1, bg2, bg3,bg4 = self.ed["Sonar/Beam_group1"], self.ed["Environment"], self.ed["Vendor_specific"],self.ed["Platform"]
        
        # Paramètres physiques
        self.z_rx_e = bg3['impedance_transceiver'].values[self.channel]
        self.sampleInterval = bg1["sample_interval"].values[self.channel][0]
        self.f_0 = bg1["transmit_frequency_start"].values[self.channel][0]
        self.f_1 = bg1["transmit_frequency_stop"].values[self.channel][0]
        self.f_c = (self.f_0 + self.f_1) / 2
        self.f_n = bg1["frequency_nominal"].values[self.channel]
        self.tau = bg1["transmit_duration_nominal"].values[self.channel][0]
        self.slope = bg1["slope"].values[self.channel][0]
        self.p_tx_e = bg1["transmit_power"].values[self.channel][0]
        self.c = bg2["sound_speed_indicative"].values
        self.temperature = bg2["temperature"].values
        self.salinity = bg2["salinity"].values
        self.acidity = bg2["acidity"].values
        self.depth = bg2["depth"].values
        self.angle_sensitivity_alongship_fn = bg1["angle_sensitivity_alongship"].values[self.channel]
        self.angle_sensitivity_athwartship_fn = bg1["angle_sensitivity_athwartship"].values[self.channel]

        self.f_m = np.linspace(self.f_0, self.f_1, self.n_f_points)
        self.lambda_f_c = self.c/self.f_c
        self.lambda_m = self.c/self.f_m
        
        self.nb_pings = bg1["backscatter_r"].values.shape[1]
        
        # Filtres WBT/PC
        h_pc = self._trim_nans(bg3["PC_filter_r"].values[self.channel], 0) + 1j*self._trim_nans(bg3["PC_filter_i"].values[self.channel], 0)
        h_wbt = self._trim_nans(bg3["WBT_filter_r"].values[self.channel], 0) + 1j*self._trim_nans(bg3["WBT_filter_i"].values[self.channel], 0)
        self.filter_v = [{"h_fl_i": h_wbt, "D": bg3["WBT_decimation"].values[self.channel]},
                        {"h_fl_i": h_pc, "D": bg3["PC_decimation"].values[self.channel]}]
        
        self.nb_pings = bg1["backscatter_r"].values.shape[1]
        self.N_u = bg1["backscatter_r"].values.shape[3]

        self.y_rx_nu = (bg1["backscatter_r"].values + 1j * bg1["backscatter_i"].values)



        ping_times = bg1.ping_time
        bg4_clean = bg4.drop_duplicates(dim="time1")
        lat_interp = bg4_clean.latitude.interp(time1=ping_times).ffill("ping_time").bfill("ping_time").values
        lon_interp = bg4_clean.longitude.interp(time1=ping_times).ffill("ping_time").bfill("ping_time").values
        lat0, lon0 = lat_interp[0], lon_interp[0]
        transformer = Transformer.from_crs("epsg:4326", "epsg:2154", always_xy=True)
        self.x_boat_gps, self.y_boat_gps = transformer.transform(lon_interp, lat_interp)



    def chirp(self,t): 
        a = np.pi * (self.f_1 - self.f_0) / self.tau
        b = 2 * np.pi * self.f_0
        return np.cos(a * t * t + b * t)

    def hann(self,L):
        n = np.arange(0, L, 1)
        return 0.5 * (1.0 - np.cos(2.0 * np.pi * n / (L - 1)))

    def generateIdealWindowedTransmitSignal(self):
        nsamples = int(np.floor(self.tau * self.f_s))
        t = np.linspace(0, nsamples - 1, num=nsamples) * 1 / self.f_s
        y = self.chirp(t)
        L = int(np.round(self.tau * self.f_s * self.slope * 2.0))  
        w = self.hann(L)
        N = len(y)
        w1 = w[0 : int(len(w) / 2)]
        w2 = w[int(len(w) / 2) :]
        i0 = 0
        i1 = len(w1)
        i2 = N - len(w2)
        i3 = N
        y[i0:i1] = y[i0:i1] * w1
        y[i2:i3] = y[i2:i3] * w2
        return y

    def calcDecmiatedSamplingRate(self):
        f_s_dec = [self.f_s]        
        if self.filter_v is not None:
            for v, _filter_v in enumerate(self.filter_v):
                f_s_dec.append(f_s_dec[v] / _filter_v["D"])
        return f_s_dec[-1]

    def calcNormalizedTransmitSignal(self):
        y_tx_n = self.generateIdealWindowedTransmitSignal()
        return y_tx_n / np.max(y_tx_n)

    def calcFilteredAndDecimatedSignal(self):
        y_tilde_tx_nv = [self.calcNormalizedTransmitSignal()]
        v = 0
        if self.filter_v is not None:
            for filter_vi in self.filter_v:
                tmp = np.convolve(y_tilde_tx_nv[v], filter_vi["h_fl_i"], mode="full")[
                    0 :: filter_vi["D"]
                ]
                y_tilde_tx_nv.append(tmp)
                v += 1
        return y_tilde_tx_nv[-1]

    def calcAutoCorrelation(self):
        y_mf_n = self.calcFilteredAndDecimatedSignal()
        f_s_dec = self.calcDecmiatedSamplingRate()
        y_mf_n_conj_rev = np.conj(y_mf_n)[::-1]
        y_mf_twoNormSquared = np.linalg.norm(y_mf_n, 2) ** 2
        y_mf_n_conj_rev = y_mf_n_conj_rev
        y_mf_twoNormSquared = y_mf_twoNormSquared
        y_mf_auto_n = np.convolve(y_mf_n, y_mf_n_conj_rev) / y_mf_twoNormSquared
        p_tx_auto = np.abs(y_mf_auto_n) ** 2
        tau_eff = np.sum(p_tx_auto) / ((np.max(p_tx_auto)) * f_s_dec)
        return y_mf_auto_n, tau_eff
    
    def calcRange(self,sampleCount):
        dr = self.sampleInterval * self.c * 0.5
        r =np.arange(sampleCount) * dr
        r[r == 0] = 1e-20 # Avoid problems with log10 for r=0
        return r

    def calcAbsorption(self,f):
        f = f / 1000

        a1 = (8.86 / self.c) * 10 ** (0.78 * self.acidity - 5)
        p1 = 1
        f1 = 2.8 * (self.salinity / 35) ** 0.5 * 10 ** (4 - 1245 / (self.temperature + 273))

        a2 = 21.44 * (self.salinity / self.c) * (1 + 0.025 * self.temperature)
        p2 = 1 - 1.37e-4 * self.depth + 6.62e-9 * self.depth**2
        f2 = 8.17 * 10 ** (8 - 1990 / (self.temperature + 273)) / (1 + 0.0018 * (self.salinity - 35))

        p3 = 1 - 3.83e-5 * self.depth + 4.9e-10 * self.depth**2

        a3l = 4.937e-4 - 2.59e-5 * self.temperature + 9.11e-7 * self.temperature**2 - 1.5e-8 * self.temperature**3
        a3h = 3.964e-4 - 1.146e-5 * self.temperature + 1.45e-7 * self.temperature**2 - 6.5e-10 * self.temperature**3
        a3 = a3l * (self.temperature <= 20) + a3h * (self.temperature > 20)

        a = f**2 * (a1 * p1 * f1 / (f1**2 + f**2)+ a2 * p2 * f2 / (f2**2 + f**2)+ a3 * p3)

        return a / 1000

    def calc_gamma_alongship(self):
        return self.angle_sensitivity_alongship_fn * (self.f_c / self.f_n)

    def calc_gamma_athwartship(self):
        return self.angle_sensitivity_athwartship_fn * (self.f_c / self.f_n)

    def calc_angle_offsets(self,f):
        angle_offset_alongship = np.interp(f, self.frequencies, self.angle_offset_alongship_fn)
        angle_offset_athwartship = np.interp(f, self.frequencies, self.angle_offset_athwartship_fn)
        return angle_offset_alongship, angle_offset_athwartship

    def calc_beam_widths(self, f):
        beam_width_alongship = np.interp(f, self.frequencies, self.beam_width_alongship_fn)
        beam_width_athwartship = np.interp(f, self.frequencies, self.beam_width_athwartship_fn)
        return beam_width_alongship, beam_width_athwartship

    def calc_b_theta_phi(self,f,theta,phi):
        angle_offset_alongship_m, angle_offset_athwartship_m = self.calc_angle_offsets(f)
        beam_width_alongship_m, beam_width_athwartship_m = self.calc_beam_widths(f)
        B_theta_phi_m = (0.5 * 6.0206 * ((np.abs(theta-angle_offset_alongship_m) / (beam_width_alongship_m / 2))** 2
                + (np.abs(phi-angle_offset_athwartship_m) / (beam_width_athwartship_m / 2))** 2
                - 0.18 * ((np.abs(theta-angle_offset_alongship_m) / (beam_width_alongship_m / 2))** 2
                    * (np.abs(phi-angle_offset_athwartship_m) / (beam_width_athwartship_m / 2))** 2)))
        return B_theta_phi_m

    def calcg0(self,f):
        dB_G0 = np.interp(f,self.frequencies,self.gain)
        return dB_G0

    def calc_g(self,f,theta,phi):
        b_theta_phi_m = self.calc_b_theta_phi(f,theta,phi)
        g0_m = self.calcg0(f)
        return 10**((g0_m - b_theta_phi_m)/10)
    
    def calcPulseCompressedSignals(self,y_rx):
        y_mf_n = self.calcFilteredAndDecimatedSignal()
        y_mf_n_conj_rev = np.conj(y_mf_n)[::-1]
        y_mf_twoNormSquared = np.linalg.norm(y_mf_n) ** 2
        pulseCompressedQuadrants = []
        start_idx = len(y_mf_n_conj_rev) - 1
        for u in y_rx:
            y_pc_nu = np.convolve(y_mf_n_conj_rev, u, mode="full") / y_mf_twoNormSquared
            y_pc_nu = y_pc_nu[start_idx::]
            pulseCompressedQuadrants.append(y_pc_nu)
        return np.array(pulseCompressedQuadrants)

    def calcAverageSignal(self,y_pc_nu):
        return np.sum(y_pc_nu, axis=0) / y_pc_nu.shape[0]

    def calcTransducerHalves(self,y_pc_nu):
        y_pc_fore_n = 0.5 * (y_pc_nu[2, :] + y_pc_nu[3, :])
        y_pc_aft_n = 0.5 * (y_pc_nu[0, :] + y_pc_nu[1, :])
        y_pc_star_n = 0.5 * (y_pc_nu[0, :] + y_pc_nu[3, :])
        y_pc_port_n = 0.5 * (y_pc_nu[1, :] + y_pc_nu[2, :])

        return y_pc_fore_n, y_pc_aft_n, y_pc_star_n, y_pc_port_n

    def calcPower(self,y_pc):
        K1 = self.N_u / ((2 * np.sqrt(2)) ** 2)
        K2 = (np.abs(self.z_rx_e + self.z_td_e) / self.z_rx_e) ** 2
        K3 = 1.0 / np.abs(self.z_td_e)
        C1Prx = K1 * K2 * K3
        Prx = C1Prx * np.abs(y_pc) ** 2
        Prx[Prx == 0] = 1e-20
        return Prx

    def calcAngles(self,y_pc_halves):
        gamma_theta,gamma_phi = self.calc_gamma_alongship(),self.calc_gamma_athwartship()
        y_pc_fore_n, y_pc_aft_n, y_pc_star_n, y_pc_port_n = y_pc_halves
        y_theta_n = y_pc_fore_n * np.conj(y_pc_aft_n)
        y_phi_n = y_pc_star_n * np.conj(y_pc_port_n)
        theta_n = (np.arcsin(np.arctan2(np.imag(y_theta_n), np.real(y_theta_n)) / gamma_theta) * 180 / np.pi)
        phi_n = (np.arcsin(np.arctan2(np.imag(y_phi_n), np.real(y_phi_n)) / gamma_phi) * 180 / np.pi)
        return theta_n, phi_n

    def calcSp(self,p_rx_e_n,theta,phi,r_n):
        S_p_n = (10.0 * np.log10(p_rx_e_n) + 40.0 * np.log10(r_n) + 2.0 * self.calcAbsorption(self.f_c)* r_n - 10* np.log10((self.p_tx_e * self.lambda_f_c**2 * self.calc_g(self.f_c,0,0)**2) / (16.0 * np.pi**2)))
        return S_p_n

    def alignAuto(self, y_pc_t_n):
        y_mf_auto_n,self.tau_eff = self.calcAutoCorrelation()
        idx_peak_auto = np.argmax(np.abs(y_mf_auto_n))
        idx_peak_target = np.argmax(np.abs(y_pc_t_n))
        n_target = len(y_pc_t_n)
        idx_start = idx_peak_auto - idx_peak_target
        idx_stop = idx_start + n_target
        y_mf_auto_red_n = y_mf_auto_n[max(0, idx_start):min(len(y_mf_auto_n), idx_stop)]
        return y_mf_auto_red_n


    def calcDFTforTS(self,y_pc_t_n, y_mf_auto_red_n, f_s_dec):
        N_DFT = int(2 ** np.ceil(np.log2(self.n_f_points)))
        idxtmp = np.floor(self.f_m / f_s_dec * N_DFT).astype("int")
        idx = np.mod(idxtmp, N_DFT)
        _Y_pc_t_m = np.fft.fft(y_pc_t_n, n=N_DFT)
        Y_pc_t_m = _Y_pc_t_m[idx]
        _Y_mf_auto_red_m = np.fft.fft(y_mf_auto_red_n, n=N_DFT)
        Y_mf_auto_red_m = _Y_mf_auto_red_m[idx]
        Y_tilde_pc_t_m = Y_pc_t_m / Y_mf_auto_red_m
        return Y_pc_t_m, Y_mf_auto_red_m, Y_tilde_pc_t_m

    def calcPowerFreqTS(self, Y_tilde_pc_t_m):
        imp = (np.abs(self.z_rx_e + self.z_td_e) / np.abs(self.z_rx_e)) ** 2 / np.abs(self.z_td_e)
        P_rx_e_t_m = self.N_u * (np.abs(Y_tilde_pc_t_m) / (2 * np.sqrt(2))) ** 2 * imp
        return P_rx_e_t_m

    def calcTSf(self,P_rx_e_t_m, r_t,theta,phi):
        TS_m = (10 * np.log10(P_rx_e_t_m) + 40 * np.log10(r_t) + 2 * self.calcAbsorption(self.f_m) * r_t - 10* np.log10((self.p_tx_e * self.lambda_m**2 * self.calc_g(self.f_m,theta,phi)**2) / (16 * np.pi**2)))
        return TS_m