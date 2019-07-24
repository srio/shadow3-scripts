import numpy
from grating_tools import m2ev

def get_sigmas_radiation(photon_energy,undulator_length):
    lambdan = m2ev / photon_energy
    return 2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),1e6*0.69*numpy.sqrt(lambdan/undulator_length),

