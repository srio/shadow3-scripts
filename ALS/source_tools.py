import numpy
from grating_tools import m2ev

def get_sigmas_radiation(photon_energy,undulator_length):
    lambdan = m2ev / photon_energy
    return 2.740/4/numpy.pi*numpy.sqrt(lambdan*undulator_length),0.69*numpy.sqrt(lambdan/undulator_length),

def get_sigmas_ALSU(photon_energy,undulator_length):
    sr, srp = get_sigmas_radiation(photon_energy, undulator_length)

    sx, sz, sxp, szp = 12.1e-6, 14.7e-6, 5.7e-6, 4.7e-6

    Sx = numpy.sqrt(sx ** 2 + sr ** 2)
    Sz = numpy.sqrt(sz ** 2 + sr ** 2)
    Sxp = numpy.sqrt(sxp ** 2 + srp ** 2)
    Szp = numpy.sqrt(szp ** 2 + srp ** 2)

    return Sx,Sz,Sxp,Szp


if __name__ == "__main__":
    energy = 806.0
    L = 2.1
    sr, srp = get_sigmas_radiation(energy, L)
    print("Radiation Sizes L=%f m at E=%f eV: \nSr=%5.3f um \nSrp=%5.3f urad: " % (L, energy, 1e6 * sr, 1e6 * srp))


    Sx, Sz, Sxp, Szp = get_sigmas_ALSU(energy,L)
    print("Sizes ALS-U L=%f m at E=%f eV: \nSx=%5.3f um \nSz=%5.3f um \nSxp=%5.3f urad \nSzp=%5.3f urad"%
          (L,energy,1e6*Sx,1e6*Sz,1e6*Sxp,1e6*Szp))