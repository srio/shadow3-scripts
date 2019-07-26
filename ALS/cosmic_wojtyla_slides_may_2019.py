

#
# follows Wojdyla slides May 2019
#
import numpy
import scipy.constants as codata

from srxraylib.plot.gol import plot


from grating_tools import m2ev, solve_grating_equation, vls_coefficients_calculate, vls_coefficients_convert_to_shadow
from grating_tools import trajectories
from source_tools import get_sigmas_ALSU
from plot_scan_results import get_shadow_result


if __name__ == "__main__":
    #
    #
    # inputs
    #

    do_plot = True
    do_overplot_shadow_results = True


    grating = 2 # 1: optimized at 379.1 eV,  2: at 806 eV,  3: at 1714 eV

    m = 1 # order

    Emin = 250.0
    Emax = 2500.0
    undulator_length = 2.1

    # grating positions
    r = 25.201
    rp = 7.573
    L = r + rp

    print("Inputs:\n  Emin=%5.3f eV,Emax=%5.3f eV "%(Emin,Emax))
    print("  r=%f, rp=%f, L=%f "%(r,rp, L))

    #
    #
    #
    if grating == 1:
        file_shadow_results = "cosmic_scan_energy_grating379eV.dat"
    elif grating == 2:
        file_shadow_results = "cosmic_scan_energy_grating806eV.dat"
    elif grating == 3:
        file_shadow_results = "cosmic_scan_energy_grating1714eV.dat"

    #
    # source size
    #

    energies = numpy.linspace(Emin, Emax, 100)

    Sx,Sz,Sxp,Szp = get_sigmas_ALSU(energies, undulator_length)

    if do_plot:
        if not do_overplot_shadow_results:
            plot(energies, Sz * 1e6, ytitle="Source Vertical sigma [um]", show=True)
        else:
            x, y = get_shadow_result(file_shadow_results, "source_size")
            plot(energies, Sz*2.35*1e6, x, y,
                 ytitle="Source Vertical FWHM [um]",
                 legend=["analytical","used in shadow"],
                 show=True)




    #
    # parameters of the used grating
    #



    if grating == 1:
        Eopt = 379.1
        # for Eopt
        Sx_opt,Sz_opt,Sxp_opt,Szp_opt = get_sigmas_ALSU(Eopt, undulator_length)

        C = 1.6320
        b2 = -0.23513474102160892 * (-1.0)
        b3 = 0.026929322694834845 * (-1.0)
        b4 = -0.0037131660182550476 * (-1.0)
        k0 = 178960
    elif grating == 2:
        Eopt = 806.0
        # for Eopt
        Sx_opt,Sz_opt,Sxp_opt,Szp_opt = get_sigmas_ALSU(Eopt, undulator_length)

        C = 1.5772
        b2 = -0.24736325719129112  * (-1.0)
        b3 = 0.028064075314566773  * (-1.0)
        b4 = -0.003883149911167572 * (-1.0)
        k0 = 287440.0
    elif grating == 3:
        Eopt = 1714.4
        # for Eopt
        Sx_opt,Sz_opt,Sxp_opt,Szp_opt = get_sigmas_ALSU(Eopt, undulator_length)

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

    size_at_source_fwhm = (2.35 * Sz)
    size_at_image_fwhm = size_at_source_fwhm * rp / r / C

    # delta_lambda0_exit = k0**(-1) / numpy.abs(m) * numpy.cos(beta) * size_at_image_fwhm / rp
    # delta_lambda0_source = k0**(-1) / numpy.abs(m) * numpy.cos(alpha) * size_at_source_fwhm / r
    # delta_lambda0 = numpy.sqrt( delta_lambda0_source**2 + delta_lambda0_exit**2)

    # print("Delta_lambda: %g, %g, %g"%(delta_lambda0_source,delta_lambda0_exit,delta_lambda0))
    # print("Resolving power: %g, %g, %g"%(wavelength/delta_lambda0_source,wavelength/delta_lambda0_exit,wavelength/delta_lambda0))


    #
    # calculate trajectories
    #


    Alpha, Beta = trajectories(energies,r,rp,k0,m,b2)
    Theta = (Alpha - Beta) / 2
    Thetap = numpy.pi - Theta
    #
    #


    if do_plot:
        if not do_overplot_shadow_results:
            plot(energies,Alpha*180/numpy.pi,energies,
                -Beta*180/numpy.pi,energies,(Alpha-Beta)*0.5*180/numpy.pi,
                xtitle="Photon energy [eV]",
                ytitle="Angle [deg]",legend=["alpha","beta","theta"])
        else:
            e, a = get_shadow_result(file_shadow_results, "alpha")
            e, b = get_shadow_result(file_shadow_results, "beta")
            plot(energies,Alpha*180/numpy.pi,
                 energies,-Beta*180/numpy.pi, #energies,(Alpha-Beta)*0.5*180/numpy.pi,
                 e,a,
                 e,b,
                xtitle="Photon energy [eV]",
                ytitle="Angle [deg]",legend=["alpha","beta","alpha (shadow)","beta (shadow)"])


    if do_plot:
        if not do_overplot_shadow_results:
            plot(energies,numpy.cos(Beta)/numpy.cos(Alpha),xtitle="Photon energy [eV]",ytitle="C")
        else:
            x, y = get_shadow_result(file_shadow_results, "C")
            plot(energies, numpy.cos(Beta) / numpy.cos(Alpha),  x, y,
                 xtitle="Photon energy [eV]", ytitle="C", legend=["analytical","used in shadow"])

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

    if do_plot:
        if not do_overplot_shadow_results:
            plot(energies, Size_at_slit*1e6,
                 xtitle="Photon energy [eV]", ytitle="Size at slit [um]",
                 yrange=[0,12])
        else:
            x, y = get_shadow_result(file_shadow_results, "image_size")
            plot(energies, Size_at_slit*1e6,
                 x, y,
                 xtitle="Photon energy [eV]", ytitle="Size at slit [um]",
                 legend=["analytical","shadow"],
                 yrange=[0,12])


    #
    # resolving power
    #

    delta_lambda_source = size_at_source_fwhm * numpy.cos(Alpha) / (k0 * m * r)
    resolving_power_source = (m2ev / energies) / delta_lambda_source

    delta_lambda_image = size_at_image_fwhm * numpy.cos(Beta) / (k0 * m * rp)
    resolving_power_image = (m2ev / energies) / delta_lambda_image


    if do_plot:
        if not do_overplot_shadow_results:
            plot(energies,resolving_power_source,
                 energies, resolving_power_image,
                 # energies, (m2ev / energies) / numpy.sqrt(delta_lambda0_source**2 + delta_lambda_image**2),
                 xtitle="Photon energy [eV]",ytitle="Resolving power",
                 legend=["source","image"])
        else:
            x, y = get_shadow_result(file_shadow_results, "resolving_power")
            x, y5 = get_shadow_result(file_shadow_results, "resolving_power5")
            plot(energies,resolving_power_source,
                 energies, resolving_power_source / numpy.sqrt(2),
                 x, y5,
                 x, y,
                 linestyle = ["solid","dashed","solid","dashed"],
                 color=["blue", "blue", "red", "red"],
                 xtitle="Photon energy [eV]",ytitle="Resolving power",
                 legend=["source","source / sqrt(2)","shadow 0.5","shadow 0.1"],yrange=[0,15000])


    if do_plot:
        y0 = Szp * 2.35 * 25.201 / numpy.cos(Alpha)
        if not do_overplot_shadow_results:
            plot(energies, y0 * 1e3,
                 yrange=[0,50],
                 xtitle="Photon energy [eV]", ytitle="Footprint on Grating [mm]")
        else:
            x, y = get_shadow_result(file_shadow_results, "footprint")
            plot(energies, y0*1e3,
                x, y * 1e3,
                 yrange=[0,50],
                xtitle="Photon energy [eV]", ytitle="Footprint on Grating [mm]",
                legend=["analytical","shadow"])