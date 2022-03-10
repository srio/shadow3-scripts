
def source(Kvalue, do_plot=0):
    #
    # script to make the calculations (created by XOPPY:undulator_spectrum)
    #
    from orangecontrib.xoppy.util.xoppy_undulators import xoppy_calc_undulator_spectrum
    energy, flux, spectral_power, cumulated_power = xoppy_calc_undulator_spectrum(
        ELECTRONENERGY=6.0,
        ELECTRONENERGYSPREAD=0.001,
        ELECTRONCURRENT=0.2,
        ELECTRONBEAMSIZEH=3.01836e-05,
        ELECTRONBEAMSIZEV=3.63641e-06,
        ELECTRONBEAMDIVERGENCEH=4.36821e-06,
        ELECTRONBEAMDIVERGENCEV=1.37498e-06,
        PERIODID=0.042,
        NPERIODS=38,
        KV=Kvalue, #1.851076,
        KH=0.0,
        KPHASE=0.0,
        DISTANCE=27.0,
        GAPH=0.001,
        GAPV=0.001,
        GAPH_CENTER=0.0,
        GAPV_CENTER=0.0,
        PHOTONENERGYMIN=2000.0,
        # PHOTONENERGYMAX=120000.0,
        # PHOTONENERGYPOINTS=5000,
        PHOTONENERGYMAX=5000.0,
        PHOTONENERGYPOINTS=500,
        METHOD=2,
        USEEMITTANCES=1)
    # example plot
    if do_plot:
        from srxraylib.plot.gol import plot
        plot(energy, flux, ytitle="Flux [photons/s/o.1%bw]", xtitle="Poton energy [eV]", title="Undulator Flux",
             xlog=False, ylog=False, show=False)
        plot(energy, spectral_power, ytitle="Power [W/eV]", xtitle="Poton energy [eV]", title="Undulator Spectral Power",
             xlog=False, ylog=False, show=False)
        plot(energy, cumulated_power, ytitle="Cumulated Power [W]", xtitle="Poton energy [eV]",
             title="Undulator Cumulated Power",
             xlog=False, ylog=False, show=True)
    #
    # end script
    #
    return energy, flux, spectral_power, cumulated_power


def optical_system(energy, spectral_power):

    # from xoppy.orangecontrib.xoppy.util.xoppy_xraylib_util import xpower_calc
    from orangecontrib.xoppy.util.xoppy_xraylib_util import xpower_calc


    xpower_calc_data =  xpower_calc(energies=energy, source=spectral_power,
                    substance=["Ni","Ni"], flags=[1,1], dens=["?","?"], thick=[0.5,0.5], angle=[6.7,6.7],
                    roughness=[0.0,0.0],output_file=None)
    return xpower_calc_data


def get_K_from_photon_energy(photon_energy,harmonic=1,period_length=0.042):
    import scipy.constants as codata
    energy_in_GeV = 6.0
    gamma = 1e9 * energy_in_GeV / (codata.m_e * codata.c ** 2 / codata.e)
    print("gamma for EBS: ", gamma)
    m2ev = codata.c * codata.h / codata.e
    wavelength = harmonic * m2ev / photon_energy
    return numpy.sqrt(2*(((wavelength*2*gamma**2)/period_length)-1))
#
# main
#

import numpy

ENERGIES = numpy.linspace(3000,15000,10)
POWER = []
K = []
HARMONIC = []

for energy in ENERGIES:

    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> energy: ", energy)
    Kvalue = get_K_from_photon_energy(energy, harmonic=1)
    if numpy.isnan(Kvalue):
        Kvalue = get_K_from_photon_energy(energy, harmonic=3)
        HARMONIC.append(3)
    else:
        HARMONIC.append(1)

    K.append(Kvalue)
    print("K @ %g eV: %g " % (energy,Kvalue))



    energy, flux, spectral_power, cumulated_power = source(Kvalue)
    xpower_calc_data = optical_system(energy, spectral_power)


    # print(xpower_calc_data["info"])
    # print(xpower_calc_data["data"].shape)

    transmitted_array = xpower_calc_data["data"][-1,:]


    power = numpy.trapz(transmitted_array, energy)

    print("Outcoming power: ", power)

    POWER.append(power)

    from srxraylib.plot.gol import plot

plot(ENERGIES, numpy.array(POWER), ytitle="Power [W]", xtitle="Poton energy [eV]", )

print("========================")
print("ENERGIES   HARMONIC   K   POWER ")
for i in range(ENERGIES.size):
    print("%d %d %g %g " % (ENERGIES[i], HARMONIC[i], K[i], POWER[i]))


