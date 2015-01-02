import scipy.constants
import pandas as pd
import numpy as np

#  ___  _____
# / _ \| ____|
#| | | |  _|
#| |_| | |___
# \__\_\_____|
#

qe_R9880U_110_par = [ 1489.053765000671, -20.61340505642701, 
        0.09607362916193821, -0.000144918944048782, -1.087924475686453e-07,
        3.619104979507752e-10, 2.742092765095943e-13, -1.067200613381487e-15,
        6.333980140159196e-19, 4.675391577876988, 505.1903283978535,
        15.37334879108591, -23.08738129086531, 358.7521218115685, 
        53.63424346389683 ]

def qe_R9880U_110( wavelength):
    par = qe_R9880U_110_par

    x = float(wavelength)
    x1 = x if x < 650 else 650
    qe = 0.0

    for i in np.arange( 0, 9 ):
        qe = qe + 1.0*par[i]*x1**float(i)

    qe =  qe + par[9]*np.exp( - 0.5 * ( (x1 - par[10] ) / par[11]) **2 )
    qe =  qe + par[12]*np.exp( - 0.5 * ( (x1- par[13] ) / par[14]) **2 )

    qe *= 0.01

    if (x>650):
        qe *= (1 - (x-650)/(675-650));

    if (qe<0 or x<200):
        qe = 0

    return float(qe)

qe_R7400U_03_par = [ -58.41145755814, 1.450540667766, -0.01561331198442,
        9.545010080831e-05, -3.648461145542e-07, 9.047599515597e-10,
        -1.457151808585e-12, 1.471328774241e-15, -8.46121819724e-19,
        2.11384701372e-22 ]

def qe_R7400U_03 (wavelength):
    par = qe_R7400U_03_par
    x = float(wavelength)

    if x<180: return float(0.0)
    if x>660: return float(0.0)

    qe = 0.0

    for i in np.arange(0, 10 ):
        qe = qe+ 1.0*par[i]*x**float(i)

    if qe<0: qe = 0.0

    return float(qe)

################################################################################

# ____ _____  _    _   _ ____    _    ____  ____  ____
#/ ___|_   _|/ \  | \ | |  _ \  / \  |  _ \|  _ \/ ___|
#\___ \ | | / _ \ |  \| | | | |/ _ \ | |_) | | | \___ \
        # ___) || |/ ___ \| |\  | |_| / ___ \|  _ <| |_| |___) |
#|____/ |_/_/   \_\_| \_|____/_/   \_\_| \_\____/|____/
#

def get_standard_eye_data():
    standard_eye = pd.read_csv( 'params/cie_v.csv', names = ['wavelength', 'eff' ] )
    return standard_eye

def get_corning_blue_data():
    corning_blue = pd.read_csv( 'params/clean_corning.dat', delim_whitespace = True )
    return corning_blue

################################################################################

# __  __    _  _____ _____ ____  ___    _    _     ____
#|  \/  |  / \|_   _| ____|  _ \|_ _|  / \  | |   / ___|
#| |\/| | / _ \ | | |  _| | |_) || |  / _ \ | |   \___ \
        #| |  | |/ ___ \| | | |___|  _ < | | / ___ \| |___ ___) |
#|_|  |_/_/   \_\_| |_____|_| \_\___/_/   \_\_____|____/
#

stp_pressure = 1.01325
stp_temperature = 273.15

def ref_n2( wavelength , pressure, temperature ):
    x = float( wavelength )
    reduced_index = 68.5520e-6 + 32431.57e-6*x*x/(144.0*x*x-1)
    pressure_factor =  pressure / stp_pressure
    temperature_factor =  temperature / stp_temperature
    reduced_index = reduced_index * pressure_factor / temperature_factor

    return 1 + reduced_index

################################################################################

# ____           _ _       _   _
#|  _ \ __ _  __| (_) __ _| |_(_) ___  _ __
#| |_) / _` |/ _` | |/ _` | __| |/ _ \| '_ \
        #|  _ < (_| | (_| | | (_| | |_| | (_) | | | |
#|_| \_\__,_|\__,_|_|\__,_|\__|_|\___/|_| |_|
#

def black_body_generator( T ):
    def bb( l ):
        h = scipy.constants.h
        c = scipy.constants.c
        k = scipy.constants.k
        print( "h {0}, c {1}, k {2}".format(  h, c, k ))
        return 2*h * c**2 / l**5 * 1 / ( np.exp(  h / l * c * 1 / ( k * T) ) - 1 )
    return bb

def cherenkov_radiation( l , b, pressure, temperature):
    return 2 * scipy.constants.pi * scipy.constants.alpha /l**2 * ( 1 - 1/b**2/ref_n2(l, pressure, temperature ) )

################################################################################

# _____                              _ _   _
#|_   _| __ __ _ _ __  ___ _ __ ___ (_) |_| |_ __ _ _ __   ___ ___
#  | || '__/ _` | '_ \/ __| '_ ` _ \| | __| __/ _` | '_ \ / __/ _ \
        #  | || | | (_| | | | \__ \ | | | | | | |_| || (_| | | | | (_|  __/
#  |_||_|  \__,_|_| |_|___/_| |_| |_|_|\__|\__\__,_|_| |_|\___\___|
#
#

def quartz_window_transmittance(iWindow, wavelength):
    x = float(wavelength)
    f = float(-999.0)

    if (iWindow==1) or (iWindow==3) or (iWindow==7):
        if (x>370):   f = 0.900 + 0.00011    * (x-370)
        elif (x>270): f = 0.910 - 0.010/100. * (x-270)
        elif (x>260): f = 0.890 + 0.020/ 10. * (x-260)
        elif (x>250): f = 0.830 + 0.060/ 10. * (x-250)
        elif (x>240): f = 0.790 + 0.040/ 10. * (x-240)
        elif (x>230): f = 0.790
        elif (x>225): f = 0.740 + 0.050/  5. * (x-225)
        elif (x>220): f = 0.640 + 0.100/  5. * (x-220)
        elif (x>210): f = 0.350 + 0.290/ 10. * (x-210)
        elif (x>205): f = 0.260 + 0.090/  5. * (x-205)
        else:         f = 0.235 + 0.025/  5. * (x-200)

    elif (iWindow==2) or (iWindow==4) or (iWindow==8):
        f = 0.9 + 0.00008*(x-200)

    elif (iWindow==5):
        if (x>380):      f = 0.92 + 0.0001*(x-380);
        elif (x>360):    f = 0.92;
        elif (x>250):    f = 0.95 - 0.030 / 110. * (x-250);
        elif (x>240):    f = 0.95;
        elif (x>210):    f = 0.925 + 0.025 / 30. * (x-210);
        else :           f = 0.900 + 0.025 / 10. * (x-200);

    elif (iWindow==6):
        if      (x>500): f = 0.946;
        elif (x>300):    f = 0.950 - 0.00002*(x-300);
        elif (x>240):    f = 0.958 - 0.008 / 60. * (x-240);
        elif (x>210):    f = 0.940 + 0.018 / 30. * (x-210);
        else:            f = 0.914 + 0.026 / 10. * (x-200);

    if (f<0.001): f = 0.001
    if (f>0.999): f = 0.999;
    return f

def mangin_mirror_reflectivity( wavelength ):
    x1   = float(wavelength)
    x  = x1 if (x1<630) else 630;
    refl = ( -0.261299
             +0.00938637*x
             -2.47446e-05*x*x
             +2.31118e-08*x*x*x
             -5.79119e-12*x*x*x*x )

    if (refl<0): refl = 0
    if (refl>1): refl = 1
    return refl

def spherical_mirror_reflectivity( wavelength ):
    x    = float(wavelength)
    refl = 0.90 if (x>450) else 0.90 - 0.05/(450-200)*(450-x)
    return refl

def external_lens_transmittance( wavelength ):
    x   = float(wavelength)
    if (x<180):  return float(0)
    if (x>=180 and x<225): return 0.500+0.01800*(x-200)
    if (x>=225 and x<235): return 0.950+0.00450*(x-225)
    if (x>=235 and x<420): return 0.995
    if (x>=420 and x<450): return 0.995-0.00050*(x-420)
    if (x>=450 and x<650): return 0.980-0.00075*(x-450)
    if (x>=650 and x<750): return 0.830-0.00030*(x-650)
    if (x>=750):          return 0.80
    return float(0)
