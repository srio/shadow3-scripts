import numpy
from srxraylib.plot.gol import plot
import scipy.constants as codata
import xraylib



def create_spectrum():
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
        PERIODID=0.018,
        NPERIODS=222,
        KV=1.76,
        KH=0.0,
        KPHASE=0.0,
        DISTANCE=27.5,
        GAPH=0.0008,
        GAPV=0.0008,
        GAPH_CENTER=0.0,
        GAPV_CENTER=0.0,
        PHOTONENERGYMIN=3000.0,
        PHOTONENERGYMAX=150000.0,
        PHOTONENERGYPOINTS=2800,
        METHOD=2,
        USEEMITTANCES=1)
    # example plot
    from srxraylib.plot.gol import plot
    plot(energy,flux,ytitle="Flux [photons/s/o.1%bw]",xtitle="Poton energy [eV]",title="Undulator Flux",
        xlog=False,ylog=False,show=False)
    plot(energy,spectral_power,ytitle="Power [W/eV]",xtitle="Poton energy [eV]",title="Undulator Spectral Power",
        xlog=False,ylog=False,show=False)
    plot(energy,cumulated_power,ytitle="Cumulated Power [W]",xtitle="Poton energy [eV]",title="Undulator Cumulated Power",
        xlog=False,ylog=False,show=True)
    #
    # end script
    #
    return energy, flux

def nist_be():
    return numpy.array([
        [1.00000E-03, 6.041E+02, 6.035E+02],
        [1.50000E-03, 1.797E+02, 1.791E+02],
        [2.00000E-03, 7.469E+01, 7.422E+01],
        [3.00000E-03, 2.127E+01, 2.090E+01],
        [4.00000E-03, 8.685E+00, 8.367E+00],
        [5.00000E-03, 4.369E+00, 4.081E+00],
        [6.00000E-03, 2.527E+00, 2.260E+00],
        [8.00000E-03, 1.124E+00, 8.839E-01],
        [1.00000E-02, 6.466E-01, 4.255E-01],
        [1.50000E-02, 3.070E-01, 1.143E-01],
        [2.00000E-02, 2.251E-01, 4.780E-02],
        [3.00000E-02, 1.792E-01, 1.898E-02],
        [4.00000E-02, 1.640E-01, 1.438E-02],
        [5.00000E-02, 1.554E-01, 1.401E-02],
        [6.00000E-02, 1.493E-01, 1.468E-02],
        [8.00000E-02, 1.401E-01, 1.658E-02],
        [1.00000E-01, 1.328E-01, 1.836E-02],
        [1.50000E-01, 1.190E-01, 2.157E-02],
        [2.00000E-01, 1.089E-01, 2.353E-02],
        [3.00000E-01, 9.463E-02, 2.548E-02],
        [4.00000E-01, 8.471E-02, 2.620E-02],
        [5.00000E-01, 7.739E-02, 2.639E-02],
        [6.00000E-01, 7.155E-02, 2.627E-02],
        [8.00000E-01, 6.286E-02, 2.565E-02],
        [1.00000E+00, 5.652E-02, 2.483E-02],
        [1.25000E+00, 5.054E-02, 2.373E-02],
        [1.50000E+00, 4.597E-02, 2.268E-02],
        [2.00000E+00, 3.938E-02, 2.083E-02],
        [3.00000E+00, 3.138E-02, 1.806E-02],
        [4.00000E+00, 2.664E-02, 1.617E-02],
        [5.00000E+00, 2.347E-02, 1.479E-02],
        [6.00000E+00, 2.121E-02, 1.377E-02],
        [8.00000E+00, 1.819E-02, 1.233E-02],
        [1.00000E+01, 1.627E-02, 1.138E-02],
        [1.50000E+01, 1.361E-02, 1.001E-02],
        [2.00000E+01, 1.227E-02, 9.294E-03]])

def diamond_filter(energy, flux, diamond_thickness_in_mm = 0.3):
    XRL_MU = numpy.zeros_like(energy)
    for i in range(energy.size):
        XRL_MU[i] = 3.51 * xraylib.CS_Total(xraylib.SymbolToAtomicNumber("C"), 1e-3*energy[i])
    return energy, flux * numpy.exp(- XRL_MU * diamond_thickness_in_mm * 1e-1)

def remove_points_for_pescao(x, y, ratio=1000.0):
    ymax = y.max()
    igood = numpy.argwhere(y > ymax / ratio)
    return x[igood].copy(), y[igood].copy()

if __name__ == "__main__":
    do_calculate_spectrum = True
    diamond_thickness_in_mm = 0.8


    outfile = "spectrumE.dat"
    rho = 1.848

    if do_calculate_spectrum:
        energy, flux = create_spectrum()
        energy, flux = diamond_filter(energy, flux, diamond_thickness_in_mm=diamond_thickness_in_mm)
        f = open(outfile, "w")
        for i in range(energy.size):
            f.write("%g  %g\n" % (energy[i], flux[i]))
        f.close()
        print("File %s written to disk." % outfile)

        energy_for_pescao, flux_for_pescao = remove_points_for_pescao(energy, flux)
        f = open("spectrumEF.dat", "w")
        for i in range(energy_for_pescao.size):
            f.write("%g  %g\n" % (energy_for_pescao[i], flux_for_pescao[i]))
        f.close()
        print("File %s written to disk." % "spectrumEF.dat")

    else: # just read file with spectrum
        a = numpy.loadtxt(outfile)
        energy = a[:,0]
        flux = a[:,1]


    spectral_power = flux * 1e3 * codata.e
    estep = (energy[1] - energy[0])
    integrated_power = (spectral_power.sum() * estep)
    print("integrated power", integrated_power)
    print("volumetric power", integrated_power / (0.8**2))


    #
    # NIST data
    #

    nist = nist_be()
    print(nist.shape)
    nist_interpolated = 10 ** numpy.interp(numpy.log10(energy), numpy.log10(1e6 * nist[:,0]), numpy.log10(rho * nist[:,2]))
    # plot(1e6 * nist[:, 0], nist[:, 1],
    #      1e6 * nist[:, 0], nist[:, 2],
    #      energy, nist_interpolated/rho, xlog=1, ylog=1,
    #      xtitle="Photon energy [eV]", ytitle="[cm2/g]")

    #
    # xraylib data
    #
    XRL_MU = numpy.zeros_like(energy)
    XRL_MU_E = numpy.zeros_like(energy)


    for i in range(energy.size):
        XRL_MU[i] = rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), 1e-3*energy[i])
        XRL_MU_E[i] = rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), 1e-3*energy[i])

    plot(
        1e-3 * energy, XRL_MU,
        1e-3 * energy, XRL_MU_E,
        1e-3 * energy, nist_interpolated,
        xlog=0, ylog=1, legend=["mu","mu_e","nist_e"],
        xtitle="Photon energy [keV]", ytitle="mu [cm^-1]")

    #
    # loop on thicknesses
    #
    THICKNESS_MM = numpy.concatenate( (numpy.linspace(0,1,100),numpy.linspace(1,10,50)))

    VOLUMETRIC_ABSORBED_POWER = numpy.zeros_like(THICKNESS_MM)
    VOLUMETRIC_ABSORBED_POWER_E = numpy.zeros_like(THICKNESS_MM)
    VOLUMETRIC_ABSORBED_POWER_NIST = numpy.zeros_like(THICKNESS_MM)

    for i, thickness_mm in enumerate(THICKNESS_MM):
        thickness_mm = THICKNESS_MM[i]
        absorbed_fraction = 1.0 - numpy.exp(-XRL_MU * thickness_mm * 1e-1)
        absorbed_fraction_e = 1.0 - numpy.exp(-XRL_MU_E * thickness_mm * 1e-1)
        absorbed_fraction_nist = 1.0 - numpy.exp(-nist_interpolated * thickness_mm * 1e-1)

        # plot(energy, absorbed_fraction, energy, absorbed_fraction_e)

        absorbed_power = (flux * absorbed_fraction * codata.e * 1e3).sum() * estep
        volumetric_absorbed_power = absorbed_power / (0.8 * 0.8 * thickness_mm)

        absorbed_power_e = (flux * absorbed_fraction_e * codata.e * 1e3).sum() * estep
        volumetric_absorbed_power_e = absorbed_power_e / (0.8 * 0.8 * thickness_mm)

        absorbed_power_nist = (flux * absorbed_fraction_nist * codata.e * 1e3).sum() * estep
        volumetric_absorbed_power_nist = absorbed_power_nist / (0.8 * 0.8 * thickness_mm)

        VOLUMETRIC_ABSORBED_POWER[i] = volumetric_absorbed_power
        VOLUMETRIC_ABSORBED_POWER_E[i] = volumetric_absorbed_power_e
        VOLUMETRIC_ABSORBED_POWER_NIST[i] = volumetric_absorbed_power_nist

        print(integrated_power, absorbed_power, volumetric_absorbed_power)
        print(integrated_power, absorbed_power_e, volumetric_absorbed_power_e)


    #
    # load pescao results and make final plot
    #
    pescao = numpy.loadtxt("pescao_0p8.dat", skiprows=2)

    plot(THICKNESS_MM, VOLUMETRIC_ABSORBED_POWER,
         THICKNESS_MM, VOLUMETRIC_ABSORBED_POWER_E,
         THICKNESS_MM, VOLUMETRIC_ABSORBED_POWER_NIST,
         pescao[:,0], pescao[:,1]/(pescao[:,0] * 0.8 * 0.8),
         xtitle="Depth [mm]", ytitle="Volumetric absorption [W/mm3]",
         title="diamond window thickness = %g mm" % diamond_thickness_in_mm,
         legend=["mu","mu_e","nist_e","Monte Carlo"])