#
# example of using xraylib to get crystal data
#
# see: 
#       http://dx.doi.org/10.1016/j.sab.2011.09.011   (paper) 
#       https://github.com/tschoonj/xraylib/          (code) 
#       http://lvserver.ugent.be/xraylib-web          (web interface, but crystals not included!)
#



#
#  import block 
#
import xraylib
import numpy as np
import scipy
import scipy.constants as codata
from scipy.special import jv

def get_wavelength(ener):
    ev2meter = codata.h*codata.c/codata.e
    return ev2meter/(ener*1e3)

def get_wavevector(ener):
    return 2 * np.pi / get_wavelength(ener)

def get_chi(ener=17.0,hh=1,kk=1,ll=1,crystal="Si"):
    #
    # get crystal data for silicon crystal
    #
    cryst = xraylib.Crystal_GetCrystal(crystal)

    #
    # define miller indices and compute dSpacing
    #
    debyeWaller = 1.0
    rel_angle = 1.0  # ratio of (incident angle)/(bragg angle) -> we work at Bragg angle
    #
    # get the structure factor (at a given energy)
    #

    fH = xraylib.Crystal_F_H_StructureFactor(cryst, ener, hh, kk, ll, debyeWaller, rel_angle)

    #
    # convert structure factor in chi (or psi) = - classical_e_radius wavelength^2 fH /(pi volume)
    #
    electron_r, tmp1, tmp2 = scipy.constants.codata.physical_constants["classical electron radius"]

    # ev2meter = codata.h*codata.c/codata.e
    # wavelength = ev2meter/(ener*1e3)
    wavelength = get_wavelength(ener)
    # print("wavelength: ",wavelength)


    volume = cryst['volume'] *1e-10*1e-10*1e-10 # volume of silicon unit cell in m^3
    cte = - electron_r * wavelength*wavelength/(np.pi * volume)
    chiH = cte*fH

    return chiH

if __name__ == "__main__":

    from srxraylib.plot.gol import plot

    print(">>>> 1 1 1,",get_chi(17.0,1,1,1,"Si"))
    print(">>>>-1-1-1,",get_chi(17.0,-1,-1,-1,"Si"))
    print(">>>> 0 0 0,",get_chi(17.0,0,0,0,"Si"))

    ener = 17.0
    alpha = 0.0

    chih    = get_chi(ener,1,1,1,"Si")
    chizero = get_chi(ener,0,0,0,"Si")
    chihm   = get_chi(ener,-1,-1,-1,"Si")

    t = 300e-6
    p = 39.0
    rayon = 1.750
    teta = xraylib.Bragg_angle(xraylib.Crystal_GetCrystal("Si"),ener,1,1,1)
    print("teta",teta*180.0/np.pi)

    Z = get_wavevector(ener) * np.sqrt( chih * chihm) / np.sin(2*teta)

    a = t * np.sin(2*teta) / (2 * np.cos(alpha))

    print("|Za|= %f >>? 1"%(np.abs(Z*a)))
    print("a = %f " % (a))

    eta = np.linspace(-1*a,1*a,100)
    R0 = jv(0,Z*np.sqrt(a**2-eta**2))
    print(eta,R0)


    q_dyn = get_wavevector(ener) * a / Z.real

    print("q_dyn: ",q_dyn)

    plot(eta,np.abs(R0)**2)







