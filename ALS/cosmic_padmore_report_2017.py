

#
# follows Padmore-Wojdyla report 4/17/2019
#
import numpy

from grating_tools import m2ev, solve_grating_equation, vls_coefficients_calculate, vls_coefficients_convert_to_shadow
from grating_tools import trajectories


if __name__ == "__main__":
    #
    # inputs
    #

    do_plot = True

    m = 1 # order

    Re = 5000.
    Rs = 10000.0

    Emin = 260.0
    Emax = 2500.0

    Eopt =numpy.sqrt(Emin * Emax)

    print("Emin=%5.3f eV,Emax=%5.3f eV,Eopt=%5.3f eV: "%(Emin,Emax,Eopt))



    Sx = 17.56e-6
    Sy = 19.45e-6

    Sxp = 20.12e-6
    Syp = 19.85e-6

    size_at_source_fwhm = 2.35 * Sy


    size_at_image_fwhm = 15e-6
    thetap = 1.25 * numpy.pi / 180.0

    L = 31.589

    #
    # calculate grating positions (eq. 24)
    #

    r = (size_at_source_fwhm * Rs / thetap) * \
        ( (L * thetap - size_at_image_fwhm * Re) / (size_at_source_fwhm * Rs + size_at_image_fwhm * Re) )

    rp = (size_at_image_fwhm * Re / thetap) * \
        ( (L * thetap + size_at_source_fwhm * Rs) / (size_at_source_fwhm * Rs + size_at_image_fwhm * Re) )

    print("r=%f, rp=%f"%(r,rp))

    #
    # calculate d-spacing
    #

    wavelength = m2ev / Eopt
    print("Wavelength= %f A"%(wavelength*1e10))
    d_spacing = wavelength / 2 / thetap**2 * \
                (size_at_source_fwhm * Rs + 2 * L * thetap - size_at_image_fwhm * Re) / \
                (size_at_source_fwhm * Rs + size_at_image_fwhm * Re)

    print("d_spacing = %g, line_density= %f"%(d_spacing,1.0/d_spacing))

    #
    # Calculate C (eq 23) and phi (eq 10)
    #

    C = (size_at_source_fwhm * Rs + L * thetap) / \
        (L * thetap - size_at_image_fwhm * Re)



    phi = thetap * ( (C - 1) / (C + 1))



    theta = 0.5 * numpy.pi - thetap

    alpha = theta + phi
    beta = -1 * (theta - phi)

    print("C=%f, phi=%f deg, alpha=%f deg, beta=%f deg" % (C,
                                                           phi * 180 / numpy.pi,
                                                           alpha * 180 / numpy.pi,
                                                           beta * 180 / numpy.pi
                                                           ))

    print("\n\n\n")


    #
    # magnification and mirror and grating positions
    #

    Mm = size_at_image_fwhm / size_at_source_fwhm


    Mg = rp / r

    print("\n\n\nL=%f m, Mm=%f"%(L,Mm))
    print("G p=%f m, q=%f m, Mg=%f"%(r,rp,Mg))


    #
    # get angles
    #

    print("\nalpha=%f deg, beta=%f deg, theta=%f deg"%(
        alpha*180/numpy.pi,beta*180/numpy.pi,theta*180/numpy.pi))


    k0 = 1.0 / d_spacing
    #
    # get resolving power
    #

    delta_lambda0_exit = k0**(-1) / numpy.abs(m) * numpy.cos(beta) * size_at_image_fwhm / rp
    delta_lambda0_source = k0**(-1) / numpy.abs(m) * numpy.cos(alpha) * (2.35 * Sy) / r
    delta_lambda0 = numpy.sqrt( delta_lambda0_source**2 + delta_lambda0_exit**2)

    print("Delta_lambda: %g, %g, %g"%(delta_lambda0_source,delta_lambda0_exit,delta_lambda0))
    print("Resolving power: %g, %g, %g"%(wavelength/delta_lambda0_source,wavelength/delta_lambda0_exit,wavelength/delta_lambda0))


    #
    # calculate grating coefficients:
    #

    b2,b3,b4 = vls_coefficients_calculate(numpy.sin(alpha), numpy.sin(beta), r, rp,
                                          line_density=k0, wavelength=wavelength, order=m)


    print(" b2=%f\n b3=%f\n b4=%f\n"%(b2,b3,b4))


    sh0,sh1,sh2,sh3 = vls_coefficients_convert_to_shadow(k0,b2,b3,b4)

    print("VLS coefficients (SI):\n k0=%f\n b2=%f\n b3=%f\n b4=%f\n"%(k0,b2,b3,b4))
    print("VLS Shadow coefficients (SI):\n c1=%f\n c2=%f\n c3=%f\n c4=%f\n" % (sh0,sh1,sh2,sh3))


    energies = numpy.linspace(Emin,Emax,100)
    Alpha, Beta = trajectories(energies,r,rp,k0,m,b2)
    Theta = (Alpha - Beta) / 2
    Thetap = numpy.pi - Theta


    from srxraylib.plot.gol import plot

    if do_plot:
        plot(energies,Alpha*180/numpy.pi,energies,-Beta*180/numpy.pi,energies,(Alpha-Beta)*0.5*180/numpy.pi,
         xtitle="Photon energy [eV]",ytitle="Angle [deg]",legend=["alpha","beta","theta"])

    if do_plot:
        plot(energies,numpy.cos(Beta)/numpy.cos(Alpha),xtitle="Photon energy [eV]",ytitle="C")

    #
    # calculate coma
    #
    w = 2 * r * numpy.sqrt(wavelength / 2 / 4.0) / numpy.cos(alpha)
    coma = 3 * w**2 * rp / ( 2 * numpy.cos(beta)) * \
           (numpy.sin(Alpha) * numpy.cos(Alpha)**2 / r**2 + \
            numpy.sin(Beta) * numpy.cos(Beta)**2 / rp**2 - \
            2 * b3 * m * k0 * m2ev / energies)

    if do_plot:
        plot(energies,coma,yrange=[-1e-8,4e-8],xrange=[200,1400],xtitle="Photon energy [eV]",ytitle="Ray deviation [m]")


    #
    # Reflectivity
    #

    from orangecontrib.xoppy.util.xoppy_xraylib_util import f1f2_calc

    reflectivity = energies * 0.0

    for i,energy in enumerate(energies):

        a = f1f2_calc("Au",energy,theta=numpy.pi/2-Alpha[i],F=8,density=None,rough=0.0,)
        reflectivity[i] = a
        # print(">>",energy,90.0-Alpha[i]*180.0/numpy.pi,a)

    if do_plot:
        plot(energies,reflectivity,yrange=[0,1],xtitle="Photon energy [eV]",ytitle="Reflectivity(Pi-Alpha(E))")

    #
    # resolving power
    #

    resolving_power_source = size_at_source_fwhm * numpy.cos(Alpha) / (k0 * m * r)
    resolving_power_source = (m2ev / energies) / resolving_power_source

    if do_plot:
        plot(energies,resolving_power_source,xtitle="Photon energy [eV]",ytitle="Resolving power (source)")