#PSP_boundaryLayerConductance.py
from __future__ import division
import numpy as np

g = 9.81
zeroKelvin = 273.15

def Kh(u, T, Tk0, atmPressure):
    cp = 29.3
    h = 0.01
    d = 0.77 * h
    zm = 0.13 * h
    zh = 0.2 * zm
    z = 1.5
    Psi_m = 0; Psi_h = 0
    vk = 0.4
    # molar density of the gas
    ro = 44.6 * (atmPressure / 101.3) * (293.15 / Tk0) 
    # volumetric heat of air (1200 J/m^3*K at 20C e sea level)
    Ch = ro * cp 
    for i in range(3):
        ustar = vk * u / (np.log((z - d + zm) / zm) + Psi_m)
        K = vk * ustar / (np.log((z - d + zh) / zh) + Psi_h)
        H = K * (T - (Tk0 - zeroKelvin))
        zeta = -vk * z * g * H / (Ch * Tk0 * np.power(ustar,3))
        if (zeta > 0):
            Psi_h = 4.7 * zeta
            Psi_m = Psi_h
        else:
            Psi_h = -2 * np.log((1 + np.sqrt(1 - 16 * zeta)) / 2)
            Psi_m = 0.6 * Psi_h
    return K
