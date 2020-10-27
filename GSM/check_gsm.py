from syned.storage_ring.light_source import LightSource
from syned.storage_ring.electron_beam import ElectronBeam
from syned.storage_ring.magnetic_structures.undulator import Undulator
from syned.util.json_tools import load_from_json_file
import numpy
import scipy.constants as codata

def mu_vartanyants(sigma,sigmaprime,wavelength):
    epsilon = sigmaprime * sigma
    kwave = 2 * numpy.pi / wavelength
    tmp = 4 * kwave**2 * epsilon**2 - 1
    return 2 * sigma / numpy.sqrt(tmp)

def beta_manolo(cf):
    # https://www.wolframalpha.com/input/?i=Solve%5B+%281-%281%2F+%281+%2B+beta**2+%2F+2%2B+beta+*+Sqrt%281%2B+%28beta%2F2%29**2%29%29++%29%29+%3Dc%2C+beta%5D
    return cf / numpy.sqrt(1-cf)

def get_coherent_fraction_exact(beta):
    q = 1 + 0.5 * beta**2 + beta * numpy.sqrt( (beta/2)**2 + 1 )
    q = 1.0 / q
    CF = 1 - q
    return CF

def get_sigmas_radiation(photon_energy,period_length=0.020, number_of_periods=100):
    # calculate sizes of the photon undulator beam
    # see formulas 25 & 30 in Elleaume (Onuki & Elleaume)
    lambdan = 1e-10 * codata.h * codata.c / codata.e * 1e10 / photon_energy  # in m
    sigma_r = 2.740 / 4 / numpy.pi * numpy.sqrt(lambdan * number_of_periods * period_length)
    sigma_r_prime =  0.69 * numpy.sqrt(lambdan / ( period_length * number_of_periods))
    return sigma_r, sigma_r_prime

def do_calculation(photon_energy=10000.0):

    a = load_from_json_file("/users/srio/Oasys/tmp_id18_u20.json")
    print(a.info())

    ebeam = a.get_electron_beam()

    mstruct = a.get_magnetic_structure()

    print(a, ebeam, mstruct)

    gamma = ebeam.gamma()
    print("gamma: ", gamma)

    K = mstruct.get_K_from_photon_energy(photon_energy, gamma,harmonic=1)
    print(">>>>>>>>>>>K: ", K)
    mstruct._K_vertical = K

    # Sx, Sz, Sxp, Szp = mstruct.get_photon_sizes_and_divergences(ebeam)
    sx, sxp, sz, szp = ebeam.get_sigmas_all()
    su, sup = get_sigmas_radiation(photon_energy,
                                   period_length=mstruct.period_length(),
                                   number_of_periods=mstruct.number_of_periods())
    Sx = numpy.sqrt( sx**2 + su**2)
    Sxp = numpy.sqrt(sxp ** 2 + sup ** 2)
    Sz = numpy.sqrt(sz ** 2 + su ** 2)
    Szp = numpy.sqrt(szp ** 2 + sup ** 2)


    cf_h = su * sup / (Sx * Sxp) # mstruct.approximated_coherent_fraction_horizontal(ebeam,harmonic=1)
    cf_v = su * sup / (Sz * Szp) # mstruct.approximated_coherent_fraction_vertical(ebeam,harmonic=1)

    wavelength = 1e-10 * codata.h * codata.c / codata.e * 1e10 / photon_energy # mstruct.resonance_wavelength(gamma)
    print("wavelength %f A" % (1e10 * wavelength))
    print("K,Sx,Sz,Sxp,Szp: ", K,Sx,Sz,Sxp,Szp )
    print("CF h,v: ", cf_h, cf_v )


    #
    # Vartanyants mu
    #
    mu_h_vartanyants = mu_vartanyants(Sx, Sxp, wavelength)
    mu_v_vartanyants = mu_vartanyants(Sz, Szp, wavelength)

    beta_h = mu_h_vartanyants / Sx
    beta_v = mu_v_vartanyants / Sz
    print("---- Vartanyants, photon energy=%f eV" %photon_energy)
    print("H sigmaI=%6.2f, sigmaMu=%6.2f, beta=%6.2f, CF_ratio=%6.2f, CF_GSM=%6.2f: " % ( 1e6*Sx, 1e6*mu_h_vartanyants, beta_h, cf_h, get_coherent_fraction_exact(beta_h)))
    print("V sigmaI=%6.2f, sigmaMu=%6.2f, beta=%6.2f, CF_ratio=%6.2f, CF_GSM=%6.2f: " % ( 1e6*Sz, 1e6*mu_v_vartanyants, beta_v, cf_v, get_coherent_fraction_exact(beta_v)))

    #
    # Manolo mu
    #
    beta_h_manolo = beta_manolo(cf_h)
    beta_v_manolo = beta_manolo(cf_v)

    mu_h_manolo = beta_h_manolo * Sx
    mu_v_manolo = beta_v_manolo * Sz

    print("---- Manolo, photon energy=%f eV" %photon_energy)
    print("H sigmaI=%6.2f, sigmaMu=%6.2f, beta=%6.2f, CF_ratio=%6.2f, CF_GSM=%6.2f: " % ( 1e6*Sx, 1e6*mu_h_manolo, beta_h_manolo, cf_h, get_coherent_fraction_exact(beta_h_manolo)))
    print("V sigmaI=%6.2f, sigmaMu=%6.2f, beta=%6.2f, CF_ratio=%6.2f, CF_GSM=%6.2f: " % ( 1e6*Sz, 1e6*mu_v_manolo, beta_v_manolo, cf_v, get_coherent_fraction_exact(beta_v_manolo)))

    return mu_h_manolo, mu_v_manolo, mu_h_vartanyants, mu_v_vartanyants

if __name__ == "__main__":
    from srxraylib.plot.gol import plot
    photon_energies = numpy.linspace(5000.0, 150000.0, 100)
    MU_H_MANOLO = numpy.zeros_like(photon_energies)
    MU_V_MANOLO = numpy.zeros_like(photon_energies)
    MU_H_VARTANYANTS = numpy.zeros_like(photon_energies)
    MU_V_VARTANYANTS = numpy.zeros_like(photon_energies)

    for i, photon_energy in enumerate(photon_energies):
        mu_h_manolo, mu_v_manolo, mu_h_vartanyants, mu_v_vartanyants = do_calculation(photon_energy=photon_energy)
        MU_H_MANOLO[i] = mu_h_manolo
        MU_V_MANOLO[i] = mu_v_manolo
        MU_H_VARTANYANTS[i] = mu_h_vartanyants
        MU_V_VARTANYANTS[i] = mu_v_vartanyants


    plot(photon_energies, 1e6*MU_H_MANOLO, photon_energies, 1e6*MU_H_VARTANYANTS, xtitle="Photon energy [eV]", ytitle="H SigmaMu",
         title="Horizontal", legend=["Manolo","Vartanyants"])
    plot(photon_energies, 1e6*MU_V_MANOLO, photon_energies, 1e6*MU_V_VARTANYANTS, xtitle="Photon energy [eV]", ytitle="V SigmaMu",
         title="Vertical", legend=["Manolo", "Vartanyants"])