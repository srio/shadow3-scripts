import numpy

from srxraylib.plot.gol import plot, set_qt
import scipy.constants as codata

set_qt()

a = numpy.loadtxt("Dia")

import xraylib

print(">>>>", a.shape)
b = numpy.zeros((a.shape[0], 3))

for i in range(a.shape[0]):
    b[i,0] = a[i,0]
    n = xraylib.Refractive_Index("C", a[i,0] * 1e-3, 3.51)
    b[i,1] = 1.0 - n.real
    b[i,2] = n.imag



plot(a[:, 0], a[:, 1],
     a[:, 0], a[:, 2],
     b[:, 0], b[:, 1],
     b[:, 0], b[:, 2],
     legend=["srcalc delta",
             "srcalc beta",
             "xraylib delta",
             "xraylib beta"],
     xlog=1,ylog=1,
     xtitle='Photon energy [eV]', ytitle="delta or beta")


# f = open("Dia-xraylib",'w')
# for i in range(a.shape[0]):
#     f.write("%g %g %g\n" % (b[i,0], b[i,1], b[i,2]) )
#
# f.close()
# print("File written to disk: Dia-xraylib")

#
#
#
density = 3.51
energy = a[:,0]
Z = 6


atwt = xraylib.AtomicWeight(Z)
avogadro = codata.Avogadro
toangstroms = codata.h * codata.c / codata.e * 1e10
re = codata.e ** 2 / codata.m_e / codata.c ** 2 / (4 * numpy.pi * codata.epsilon_0) * 1e2  # in cm




molecules_per_cc = density * avogadro / atwt
wavelength = toangstroms / energy * 1e-8  # in cm
k = molecules_per_cc * re * wavelength * wavelength / 2.0 / numpy.pi

f1 = numpy.zeros_like(energy)
f2 = numpy.zeros_like(energy)
mu = numpy.zeros_like(energy)
mu_mass = numpy.zeros_like(energy)

# elif F == 5:  # F=5  returns Photoelectric linear absorption coefficient
# out = numpy.zeros_like(energy)
# for i, ienergy in enumerate(energy):
#     out[i] = density * xraylib.CS_Photo(Z, 1e-3 * ienergy)
# elif F == 6:  # F=6  returns Photoelectric mass absorption coefficient
# out = numpy.zeros_like(energy)
# for i, ienergy in enumerate(energy):
#     out[i] = xraylib.CS_Photo(Z, 1e-3 * ienergy)

# mu = 4.0 * numpy.pi * beta / wavelength

for i, ienergy in enumerate(energy):
    f1[i] = Z + xraylib.Fi(Z, 1e-3 * ienergy)
    f2[i] = - xraylib.Fii(Z, 1e-3 * ienergy)
    mu[i] = density * xraylib.CS_Total(Z, 1e-3 * ienergy)
    mu_mass[i] =density * xraylib.CS_Energy(Z, 1e-3 * ienergy)

wavelength = codata.h * codata.c / codata.e / (energy)
mu_beta = 4.0 * numpy.pi * b[:,2] / (wavelength * 1e2)
beta_mass = mu_mass / 4 / numpy.pi  * (wavelength * 1e2)

alpha = 2.0 * k * f1
gamma = 2.0 * k * f2

plot(
     a[:, 0], a[:, 2],
     b[:, 0], b[:, 2],
     energy, beta_mass,
     legend=[
             "srcalc beta",
             "xraylib beta",
             "xraylib beta mass"],
     xlog=1,ylog=1,
     xtitle='Photon energy [eV]', ytitle="delta or beta")


#
# plot(energy, mu,
#      energy, mu_mass,
#      energy, mu_beta*1.1,
#      xlog=1,ylog=1,
#      legend=["mu","mu_mass","mu_beta*1.1"])




    # Energy            μ / ρ            μen / ρ
    # (MeV)(cm2 / g) (cm 2 /g)
nist_data_MeV = numpy.array([
    1.00000E-03,
    1.50000E-03,
    2.00000E-03,
    3.00000E-03,
    4.00000E-03,
    5.00000E-03,
    6.00000E-03,
    8.00000E-03,
    1.00000E-02,
    1.50000E-02,
    2.00000E-02,
    3.00000E-02,
    4.00000E-02,
    5.00000E-02,
    6.00000E-02,
    8.00000E-02,
    1.00000E-01,
    1.50000E-01,
    2.00000E-01,
    3.00000E-01,
    4.00000E-01,
    5.00000E-01,
    6.00000E-01,
    8.00000E-01,
    1.00000E+00,
    1.25000E+00,
    1.50000E+00,
    2.00000E+00,
    3.00000E+00,
    4.00000E+00,
    5.00000E+00,
    6.00000E+00,
    8.00000E+00,
    1.00000E+01,
    1.50000E+01,
    2.00000E+01,])

nist_data_mu = density * numpy.array([
    2.211E+03,
    7.002E+02,
    3.026E+02,
    9.033E+01,
    3.778E+01,
    1.912E+01,
    1.095E+01,
    4.576E+00,
    2.373E+00,
    8.071E-01,
    4.420E-01,
    2.562E-01,
    2.076E-01,
    1.871E-01,
    1.753E-01,
    1.610E-01,
    1.514E-01,
    1.347E-01,
    1.229E-01,
    1.066E-01,
    9.546E-02,
    8.715E-02,
    8.058E-02,
    7.076E-02,
    6.361E-02,
    5.690E-02,
    5.179E-02,
    4.442E-02,
    3.562E-02,
    3.047E-02,
    2.708E-02,
    2.469E-02,
    2.154E-02,
    1.959E-02,
    1.698E-02,
    1.575E-02,])

nist_data_mu_mass = density * numpy.array([
    2.209E+03,
    6.990E+02,
    3.016E+02,
    8.963E+01,
    3.723E+01,
    1.866E+01,
    1.054E+01,
    4.242E+00,
    2.078E+00,
    5.627E-01,
    2.238E-01,
    6.614E-02,
    3.343E-02,
    2.397E-02,
    2.098E-02,
    2.037E-02,
    2.147E-02,
    2.449E-02,
    2.655E-02,
    2.870E-02,
    2.950E-02,
    2.969E-02,
    2.956E-02,
    2.885E-02,
    2.792E-02,
    2.669E-02,
    2.551E-02,
    2.345E-02,
    2.048E-02,
    1.849E-02,
    1.710E-02,
    1.607E-02,
    1.468E-02,
    1.380E-02,
    1.258E-02,
    1.198E-02, ])


nist_data_eV = nist_data_MeV * 1e6

nist_wavelength = codata.h * codata.c / codata.e / (nist_data_eV)





plot(nist_data_eV, nist_data_mu,
     nist_data_eV, nist_data_mu_mass,
     energy, mu,
     energy, mu_mass,
     legend=["mu NIST","mu mass NIST","mu xraylib","mu mass xraylib"],
     xlog=1,ylog=1)

beta_nist_mass = nist_data_mu_mass / 4 / numpy.pi  * (nist_wavelength * 1e2)

beta_nist_mass_interpolated = 10 ** numpy.interp(numpy.log10(energy),
                                                 numpy.log10(nist_data_eV),
                                                 numpy.log10(beta_nist_mass))
# nist_interpolated = 10 ** numpy.interp(numpy.log10(energy), numpy.log10(1e6 * nist[:, 0]),
#                                        numpy.log10(rho * nist[:, 2]))


plot(energy, beta_mass,
     nist_data_eV, beta_nist_mass,
     energy, beta_nist_mass_interpolated,
     energy,a[:,2],
     xlog=1,ylog=1,legend=["beta_mass xraylib", "beta_mass nist", "beta_mass nist INTERPOLATED","beta SRCALC"])

f = open("Dia-nist",'w')
for i in range(a.shape[0]):
    if b[i,0] < 7000:
        f.write("%g %g %g\n" % (a[i, 0], a[i, 1], a[i, 2]))
    else:
        f.write("%g %g %g\n" % (a[i,0], a[i,1], beta_nist_mass_interpolated[i]) )

f.close()
print("File written to disk: Dia-nist")

c = numpy.loadtxt("Dia-nist")

# plot(
#     a[:, 0], a[:, 1],
#     c[:, 0], c[:, 1],
#     a[:, 0], a[:, 2],
#     c[:, 0], c[:, 2],
#     xlog=1,ylog=1,legend=['srcalc delta','new delta','srcalc beta','new beta'],
#     )

plot(
    a[:, 0], a[:, 2],
    c[:, 0], c[:, 2],
    b[:, 0], b[:, 2],
    xlog=1,ylog=1,legend=['srcalc beta','new beta','yesterday beta'],
    )