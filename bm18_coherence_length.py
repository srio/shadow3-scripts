import numpy as np
from srxraylib.plot.gol import plot
import scipy.constants as codata



#inputs BM18-EBS
distance_bm = 180.0 # 220.0
beta_bm = 1.814
eta_bm = 0.018

#inputs ID19-EBS
distance_und = 155.0
beta_und = 6.9
eta_und = 0.0017
undulator_length = 2.0

#inputs ID19-ESRF
distance_und_present = 155.0
beta_und_present = 0.348
eta_und_present = 0.031



emittance_ebs = 134e-12
emittance_present = 4000e-12
energy_dispersion = 0.001

photon_energy = np.linspace(20000.,270000.,200)
wavelength = codata.h * codata.c / photon_energy / codata.e

#
# BM18
#
sigma_bm = np.sqrt( emittance_ebs * beta_bm + (eta_bm * energy_dispersion)**2 )
print("sigma_bm: %f, corrected by dispersion: %f"%(1e6*np.sqrt( emittance_ebs * beta_bm),1e6*sigma_bm))
size_broadening = 0.0 # 34e-6 # transversal excursion of electrons in the short wiggler
size_bm = 2.355*sigma_bm + size_broadening # just addition, not in quadrature as size_broadening is not Gaussian
coherence_length_bm = wavelength * distance_bm / ( 2 * size_bm)

#
# ID19-EBS
#
sigma_und = np.sqrt( emittance_ebs * beta_und + (eta_und * energy_dispersion)**2 )
print("sigma_und: %f, corrected by dispersion: %f"%(1e6*np.sqrt( emittance_ebs * beta_und),1e6*sigma_und))
sigma_phot = np.sqrt(wavelength / 2 / undulator_length)
size_und = 2.355 * np.sqrt( (sigma_und)**2 + (sigma_phot)**2 )
coherence_length_und = wavelength * distance_und / ( 2 * size_und)

#
# ID19-ESRF
#
sigma_und_present = np.sqrt( emittance_present * beta_und_present + (eta_und_present * energy_dispersion)**2 )
print("sigma_und_present: %f, corrected by dispersion: %f"%(1e6*np.sqrt( emittance_present * beta_und_present),1e6*sigma_und_present))
sigma_phot = np.sqrt(wavelength / 2 / undulator_length)
size_und_present = 2.355 * np.sqrt( (sigma_und_present)**2 + (sigma_phot)**2 )
coherence_length_und_present = wavelength * distance_und / ( 2 * size_und_present)

print("size_bm: %f , size_und %f: ratio: %f, size_und_present: %f"%(1e6*size_bm,1e6*size_und[0],size_bm/size_und[0],1e6*size_und_present[0]))
print("size_bm/2.35: %f , size_und/2.35 %f: ratio: %f, size_und_present/2.35: %f"%(
     1e6*size_bm/2.35,1e6*size_und[0]/2.35,size_bm/size_und[0],1e6*size_und_present[0]/2.35))

plot(1e-3*photon_energy,1e6*coherence_length_bm,
     1e-3*photon_energy,1e6*coherence_length_und,
     1e-3*photon_energy,1e6*coherence_length_und_present,
     xlog=False,ylog=True,xtitle="photon energy [keV]",ytitle="coherence length [um]",
     yrange=[1,1000],legend=["bm18","id19-EBS","id19-present"])

# # Paul's values
#
# size_bm          = 26.0e-6  / 2.355 # BAD
# size_und         = 30.3e-6  / 2.355 # BAD
# size_und_present = 150.0e-6 / 2.355 # BAD
# size_bm          = 26.0e-6  * 2.355 # GOOD
# size_und         = 30.3e-6  * 2.355 # GOOD
# size_und_present = 150.0e-6         # GOOD
# coherence_length_bm = wavelength * distance_bm / ( 2 * size_bm)
# coherence_length_und = wavelength * distance_und / ( 2 * size_und)
# coherence_length_und_present = wavelength * distance_und / ( 2 * size_und_present)
#
# print("size_bm: %f , size_und %f: ratio: %f, size_und_present: %f"%(1e6*size_bm,1e6*size_und,size_bm/size_und,1e6*size_und_present))
# print("ration distance bm/diatance und: %f"%(distance_bm/distance_und))

