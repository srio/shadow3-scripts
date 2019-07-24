

#
# follows Wojdyla slides May 2019
#
import numpy
import scipy.constants as codata

from srxraylib.plot.gol import plot


from grating_tools import m2ev, solve_grating_equation, vls_coefficients_calculate, vls_coefficients_convert_to_shadow
from grating_tools import trajectories
from source_tools import get_sigmas_radiation


if __name__ == "__main__":
    #
    # inputs
    #

    do_plot = True

    m = 1 # order


    Emin = 250.0
    Emax = 2500.0

    print("Emin=%5.3f eV,Emax=%5.3f eV "%(Emin,Emax))


    # grating positions

    r = 25.201
    rp = 7.573
    L = r + rp
    print("r=%f, rp=%f, L=%f "%(r,rp, L))


    #
    # source size
    #

    energies = numpy.linspace(Emin, Emax, 100)
    sr, srp = get_sigmas_radiation(energies, 2.1)


    # Sx = 17.5e-6
    # Sy = 19.4e-6
    #
    # Sxp = 20.2e-6
    # Syp = 20.0e-6

    sx, sz, sxp, szp = 12.1e-6, 14.7e-6, 5.7e-6, 4.7e-6

    # these are arrays
    Sx = numpy.sqrt(sx ** 2 + sr ** 2)
    Sz = numpy.sqrt(sz ** 2 + sr ** 2)
    Sxp = numpy.sqrt(sxp ** 2 + srp ** 2)
    Szp = numpy.sqrt(szp ** 2 + srp ** 2)

    plot(energies,Sz*1e6,ytitle="Source Vertical sigma [um]")







    #
    # parameters of the used grating
    #

    grating = 3

    if grating == 1:
        Eopt = 379.1
        # for Eopt
        sr, srp = get_sigmas_radiation(Eopt, 2.1)
        Sx_opt = numpy.sqrt(sx ** 2 + sr ** 2)
        Sz_opt = numpy.sqrt(sz ** 2 + sr ** 2)
        Sxp_opt = numpy.sqrt(sxp ** 2 + srp ** 2)
        Szp_opt = numpy.sqrt(szp ** 2 + srp ** 2)

        C = 1.6320
        b2 = -0.23513474102160892 * (-1.0)
        b3 = 0.026929322694834845 * (-1.0)
        b4 = -0.0037131660182550476 * (-1.0)
        k0 = 178960
    elif grating == 2:
        Eopt = 806.0
        # for Eopt
        sr, srp = get_sigmas_radiation(Eopt, 2.1)
        Sx_opt = numpy.sqrt(sx ** 2 + sr ** 2)
        Sz_opt = numpy.sqrt(sz ** 2 + sr ** 2)
        Sxp_opt = numpy.sqrt(sxp ** 2 + srp ** 2)
        Szp_opt = numpy.sqrt(szp ** 2 + srp ** 2)

        C = 1.5772
        b2 = -0.24736325719129112  * (-1.0)
        b3 = 0.028064075314566773  * (-1.0)
        b4 = -0.003883149911167572 * (-1.0)
        k0 = 287440.0
    elif grating == 3:
        Eopt = 1714.4
        # for Eopt
        sr, srp = get_sigmas_radiation(Eopt, 2.1)
        Sx_opt = numpy.sqrt(sx ** 2 + sr ** 2)
        Sz_opt = numpy.sqrt(sz ** 2 + sr ** 2)
        Sxp_opt = numpy.sqrt(sxp ** 2 + srp ** 2)
        Szp_opt = numpy.sqrt(szp ** 2 + srp ** 2)

        C = 1.7313
        b2 = -0.21796876033899018 * (-1)
        b3 = 0.025361695845711445 * (-1)
        b4 = -0.0034823003997529007 * (-1)
        k0 = 352450.0


    #
    # get resolving power
    #

    alpha, beta = trajectories(Eopt, r, rp, k0, m, b2, verbose=True)

    print("for energy=%f eV, alpha=%f deg, beta=%f deg" % (Eopt, alpha*180/numpy.pi, beta*180/numpy.pi))
    wavelength = m2ev / Eopt

    size_at_source_fwhm = (2.35 * Sz_opt)
    size_at_image_fwhm = size_at_source_fwhm * rp / r / C

    delta_lambda0_exit = k0**(-1) / numpy.abs(m) * numpy.cos(beta) * size_at_image_fwhm / rp
    delta_lambda0_source = k0**(-1) / numpy.abs(m) * numpy.cos(alpha) * size_at_source_fwhm / r
    delta_lambda0 = numpy.sqrt( delta_lambda0_source**2 + delta_lambda0_exit**2)

    print("Delta_lambda: %g, %g, %g"%(delta_lambda0_source,delta_lambda0_exit,delta_lambda0))
    print("Resolving power: %g, %g, %g"%(wavelength/delta_lambda0_source,wavelength/delta_lambda0_exit,wavelength/delta_lambda0))


    #
    # calculate trajectories
    #


    Alpha, Beta = trajectories(energies,r,rp,k0,m,b2)
    Theta = (Alpha - Beta) / 2
    Thetap = numpy.pi - Theta
    #
    #


    if do_plot:
        plot(energies,Alpha*180/numpy.pi,energies,-Beta*180/numpy.pi,energies,(Alpha-Beta)*0.5*180/numpy.pi,
         xtitle="Photon energy [eV]",ytitle="Angle [deg]",legend=["alpha","beta","theta"])

    # if do_plot:
    #     plot(energies,numpy.cos(Beta)/numpy.cos(Alpha),xtitle="Photon energy [eV]",ytitle="C")


    #
    # #
    # # calculate coma
    # #
    # w = 2 * r * numpy.sqrt(wavelength / 2 / 4.0) / numpy.cos(alpha)
    # coma = 3 * w**2 * rp / ( 2 * numpy.cos(beta)) * \
    #        (numpy.sin(Alpha) * numpy.cos(Alpha)**2 / r**2 + \
    #         numpy.sin(Beta) * numpy.cos(Beta)**2 / rp**2 - \
    #         2 * b3 * m * k0 * m2ev / energies)
    #
    # if do_plot:
    #     plot(energies,coma,yrange=[-1e-8,4e-8],xrange=[200,1400],xtitle="Photon energy [eV]",ytitle="Ray deviation [m]")
    #
    #
    # #
    # # Reflectivity
    # #
    #
    # from orangecontrib.xoppy.util.xoppy_xraylib_util import f1f2_calc
    #
    # reflectivity = energies * 0.0
    #
    # for i,energy in enumerate(energies):
    #
    #     a = f1f2_calc("Au",energy,theta=numpy.pi/2-Alpha[i],F=8,density=None,rough=0.0,)
    #     reflectivity[i] = a
    #     # print(">>",energy,90.0-Alpha[i]*180.0/numpy.pi,a)
    #
    # if do_plot:
    #     plot(energies,reflectivity,yrange=[0,1],xtitle="Photon energy [eV]",ytitle="Reflectivity(Pi-Alpha(E))")
    #

    #
    # Beam size at slit
    #

    Size_at_slit = 2.35 * Sz * rp / r / (numpy.cos(Beta)/numpy.cos(Alpha))

    plot(energies, Size_at_slit*1e6, xtitle="Photon energy [eV]", ytitle="Size at slit [um]")

    #
    # resolving power
    #

    delta_lambda_source = size_at_source_fwhm * numpy.cos(Alpha) / (k0 * m * r)
    resolving_power_source = (m2ev / energies) / delta_lambda_source

    delta_lambda_image = size_at_image_fwhm * numpy.cos(Beta) / (k0 * m * rp)
    resolving_power_image = (m2ev / energies) / delta_lambda_image


    if do_plot:
        plot(energies,resolving_power_source,
             energies, resolving_power_image,
             energies, (m2ev / energies) / numpy.sqrt(delta_lambda0_source**2 + delta_lambda_image**2),
             xtitle="Photon energy [eV]",ytitle="Resolving power",
             legend=["source","image","combined"])