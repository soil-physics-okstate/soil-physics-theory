#PSP_soil.py

CAMPBELL = 1
VAN_GENUCHTEN = 2
waterRetentionCurve = 0

thetaS      = 0.46          #[m3/m3] saturated water content
thetaR      = 0.01          #[m3/m3] residual water content
Campbell_he = -4.2          #[J/kg] air-entry potential in Campbell formulation
Campbell_b  = 3.58          #[-] slope parameter in Campbell formulation
VG_alpha    = 0.18          #[kg/J] alpha parameter in van Genuchten formulation
VG_n        = 2.86          #[-] n parameter in van Genuchten formulation
VG_m        = 0.11          #[-] m parameter in van Genuchten formulation

def waterContent(signPsi):
    if (waterRetentionCurve == CAMPBELL):
        if (signPsi >= Campbell_he):
            return thetaS
        else:
            Se = (signPsi / Campbell_he)**(-1. / Campbell_b)
            return Se*thetaS
        
    elif(waterRetentionCurve == VAN_GENUCHTEN):
        if (signPsi >= 0.):
            return thetaS
        else:
            Se = 1. / (1. + (VG_alpha * abs(signPsi))**VG_n)**VG_m
            return Se*(thetaS - thetaR) + thetaR