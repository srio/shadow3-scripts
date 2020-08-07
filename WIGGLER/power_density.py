"""
File with trajectory written to file: /users/srio/Oasys/tmp.traj

wiggler_cdf: Electron beam energy (from velocities) = 3.000355 GeV

wiggler_cdf: gamma (from velocities) = 5870.853556 GeV
wiggler_cdf: Curvature (min)) = 0.000000 m^-1
wiggler_cdf:           (max)    0.199920 m^-1
wiggler_cdf: Radius of curvature (max) = 81689012171814624.000000 m
wiggler_cdf:                     (min) = 5.002009 m
wiggler_cdf: Critical Energy (max.) = 11973.937061 eV
wiggler_cdf:                 (min.) = 0.000000 eV
wiggler_cdf: Total no.of photons = 1.690471e+17 (in DE=99900.000 eV)
wiggler_cdf: File with wiggler cdf written to file: b'/users/srio/Oasys/xshwig.sha'

Electron beam energy (from velocities) = 3.000355 GeV

gamma (from velocities) = 5870.851896
curvature (max) = 0.199920 m
          (min) = 0.000000 m
Radius of curvature (max) = 81689012171830928.000000 m
                    (min) = 5.002009 m
Critical Energy (max.) = 11973.926903 eV
                (min.) = 0.000000 eV
File with wiggler spectrum written to file: spectrum.dat

Total power (from integral of spectrum): 10106.973910 W

Total number of photons (from integral of spectrum): 1.62115e+19

"""

#
# script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
#
from srxraylib.sources import srfunc
from srxraylib.plot.gol import plot, plot_image, plot_scatter, plot_show
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter

from srxraylib.plot.gol import set_qt
from scipy.interpolate import interp1d
set_qt()

def P(u):
    return 2 * numpy.pi / numpy.sqrt(3) * u * srfunc.fintk53(u)

def xoppy_calc_wiggler_radiation(
        ELECTRONENERGY           = 3.0,
        ELECTRONENERGYSPREAD     = 0.0,
        ELECTRONCURRENT          = 0.1,
        ELECTRONBEAMSIZEH        = 10e-6,
        ELECTRONBEAMSIZEV        = 10e-6,
        ELECTRONBEAMDIVERGENCEH  = 10e-6,
        ELECTRONBEAMDIVERGENCEV  = 10e-6,
        PERIODID                 = 0.120,
        NPERIODS                 = 37,
        KV                       = 22.416,
        KH                       = 0.0,
        KPHASE                   = 0.0,
        DISTANCE                 = 30.0,
        GAPH                     = None,
        GAPV                     = None,
        HSLITPOINTS              = 500,
        VSLITPOINTS              = 500,
        METHOD                   = 0,
        PHOTONENERGYMIN          = 100.0,
        PHOTONENERGYMAX          = 100100.0,
        PHOTONENERGYPOINTS       = 101,
        USEEMITTANCES            = 0,
        h5_file                  = "wiggler_radiation.h5",
        h5_entry_name            = "XOPPY_RADIATION",
        h5_initialize            = True,
        h5_parameters            = None,
        ):

    (traj, pars) = srfunc.wiggler_trajectory(
        b_from = 0,
        inData = "",
        nPer = NPERIODS, #37,
        nTrajPoints = HSLITPOINTS,
        ener_gev = ELECTRONENERGY,
        per = PERIODID,
        kValue = KV,
        trajFile = "tmp.traj",
        shift_x_flag = 0,
        shift_x_value = 0.0,
        shift_betax_flag = 0,
        shift_betax_value = 0.0)


    energy, flux, power = srfunc.wiggler_spectrum(traj,
        enerMin = PHOTONENERGYMIN,
        enerMax = PHOTONENERGYMAX,
        nPoints = PHOTONENERGYPOINTS,
        electronCurrent = ELECTRONCURRENT,
        outFile = "",
        elliptical = False)


    # #
    # # calculate cdf and write file for Shadow/Source
    # #
    #
    # tmp = srfunc.wiggler_cdf(traj,
    #     enerMin = 100.0,
    #     enerMax = 100000.0,
    #     enerPoints = 1001,
    #     outFile = b'tmp.sha',
    #     elliptical = False)
    #
    # print(">>>>>",tmp)

    gamma = ELECTRONENERGY / 512e-6

    X = traj[0,:].copy()
    Y = traj[1,:].copy()
    Z = traj[1,:].copy()
    divX = traj[3,:].copy()


    divZ = traj[5,:].copy()


    curX = traj[6,:].copy()
    By = traj[7, :].copy()
    # posX = divX * (distance + Y)


    Ec = 665.0 * 3**2 * numpy.abs(By)
    Ecmax = 665.0 * 3 ** 2 * numpy.abs(By.max())

    sigmaBp = 0.597 / gamma * numpy.sqrt(Ecmax / PHOTONENERGYMIN)
    divXX = numpy.linspace(divX.min() - 3 * sigmaBp, divX.max() + 3 * sigmaBp, HSLITPOINTS)
    divZZ = numpy.linspace(-3 * sigmaBp, 3 * sigmaBp, VSLITPOINTS)


    e = numpy.linspace(PHOTONENERGYMIN, PHOTONENERGYMAX, PHOTONENERGYPOINTS)

    p = numpy.zeros( (PHOTONENERGYPOINTS, HSLITPOINTS, VSLITPOINTS) )

    if PHOTONENERGYPOINTS > 3:
        do_plot = False
    else:
        do_plot = True

    for i in range(e.size):

        Ephoton = e[i]

        # horizontal divergence after Tanaka
        if False:
            e_over_ec = Ephoton / Ecmax
            uudlim = 1.0 / gamma
            print(">>>>>gamma",gamma)

            uud = numpy.linspace(-uudlim*0.99, uudlim*0.99, divX.size)
            uu  = e_over_ec / numpy.sqrt(1 - gamma**2 * uud**2)
            plot(uud, P(uu))

        # vertical divergence
        fluxDivZZ = srfunc.sync_ang(1,divZZ * 1e3,polarization=0,
               e_gev=3,i_a=0.1,hdiv_mrad=1.0,energy=Ephoton, ec_ev=Ecmax)

        if do_plot:
            plot(divZZ, fluxDivZZ, title="min intensity %f" % fluxDivZZ.min(), xtitle="divZ", ytitle="fluxDivZZ", show=1)


        # horizontal divergence
        intensity = P(Ephoton / Ec)
        fintensity = interp1d(divX, intensity, kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=0.0,
                                   assume_sorted=False)
        intensity_interpolated = fintensity(divXX)

        if True:
            intensity_interpolated.shape = -1
            fluxDivZZCC = srfunc.sync_ang(1, divXX * 1e3, polarization=0,
                                        e_gev=3, i_a=0.1, hdiv_mrad=1.0, energy=Ephoton, ec_ev=Ecmax)
            fluxDivZZCC.shape = -1

            print(">>>>>>>", intensity_interpolated.shape, fluxDivZZCC.shape)
            intensity_convolved = numpy.convolve(intensity_interpolated/intensity_interpolated.max(),
                                                 fluxDivZZCC/fluxDivZZCC.max(),
                                                 mode='same')
        else:
            intensity_convolved = intensity_interpolated


        if do_plot:
            plot(divX, intensity/intensity.max(),
                 divXX, intensity_interpolated/intensity_interpolated.max(),
                 divXX, intensity_convolved/intensity_convolved.max(),
                 title=">>>>> min intensity %f, Ephoton=%6.2f" % (intensity.min(), Ephoton), xtitle="divX", ytitle="intensity",
                 legend=["orig","interpolated","convolved"],show=1)




        # combine H * V
        INTENSITY = numpy.outer(intensity_convolved/intensity_convolved.max(), fluxDivZZ/fluxDivZZ.max())

        print(">>>>", flux.shape, INTENSITY.shape, p.shape)
        p[i,:,:] = INTENSITY / INTENSITY.sum() * flux[i]

        if do_plot:
            plot_image(INTENSITY, divXX, divZZ, aspect='auto', title="E=%6.2f" % Ephoton, show=1)
    #

    h = divXX * DISTANCE
    v = divZZ * DISTANCE


    if h5_file != "":
        try:
            if h5_initialize:
                h5w = H5SimpleWriter.initialize_file(h5_file,creator="xoppy_wigglers.py")
            else:
                h5w = H5SimpleWriter(h5_file,None)
            h5w.create_entry(h5_entry_name,nx_default=None)
            h5w.add_stack(e,h,v,p,stack_name="Radiation",entry_name=h5_entry_name,
                title_0="Photon energy [eV]",
                title_1="X gap [mm]",
                title_2="Y gap [mm]")
            h5w.create_entry("parameters",root_entry=h5_entry_name,nx_default=None)
            #TODO: open!
            # for key in h5_parameters.keys():
            #     h5w.add_key(key,h5_parameters[key], entry_name=h5_entry_name+"/parameters")
            print("File written to disk: %s"%h5_file)
        except:
            print("ERROR initializing h5 file")


    return e, h, v, p


if __name__ == "__main__":

    e, h, v, p = xoppy_calc_wiggler_radiation()


