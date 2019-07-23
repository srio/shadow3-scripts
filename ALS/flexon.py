
import scipy.constants as codata
import numpy



m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)


# # used by Luca
# codata_h = numpy.array(6.62606957e-34)
# codata_ec = numpy.array(1.602176565e-19)
# codata_c = numpy.array(299792458.0)
# m2ev = codata_c * codata_h / codata_ec      # lambda(m)  = m2eV / energy(eV)

def solve_grating_equation(line_density=100000,wavelength=10e-10,c=1.10,order=1,method='beta'):

    if method == 'beta':  # beta!!

        A = (1.0 - 1.0 / c ** 2)
        B = -2 * order * line_density * wavelength
        C = order ** 2 * line_density ** 2 * wavelength ** 2 - 1 + 1.0 / c ** 2

        Delta = B ** 2 - 4 * A * C

        sinbeta1 = (-B + numpy.sqrt(Delta)) / 2 / A
        sinbeta2 = (-B - numpy.sqrt(Delta)) / 2 / A

        # print("Discriminant=%f, sinbeta1=%f, sinbeta2=%f" % (Delta, sinbeta1, sinbeta2))

        if numpy.abs(sinbeta1) <= 1:
            sinbeta = sinbeta1
        elif numpy.abs(sinbeta2) <= 1:
            sinbeta = sinbeta2
        else:
            raise Exception("No valid beta angle")

        # change sign of beta
        # sinbeta *= -1

        # beta = numpy.arcsin(sinbeta)

        sinalpha = order * wavelength * line_density - sinbeta

        # alpha = numpy.arcsin(sinalpha)

    else:  # alpha
        A = (1.0 - order ** 2)
        B = -2 * order * line_density * wavelength
        C = order ** 2 * line_density ** 2 * wavelength ** 2 - 1 + c ** 2

        Delta = B ** 2 - 4 * A * C

        sinalpha1 = (-B + numpy.sqrt(Delta)) / 2 / A
        sinalpha2 = (-B - numpy.sqrt(Delta)) / 2 / A

        # print("Discriminant=%f, sinalpha1=%f, sinalpha2=%f" % (Delta, sinalpha1, sinalpha2))

        if numpy.abs(sinalpha1) <= 1:
            sinalpha = sinalpha1
        elif numpy.abs(sinalpha2) <= 1:
            sinalpha = sinalpha2
        else:
            raise Exception("No valid alpha angle")

        # # my value
        # sinalpha_me = m*lambda0*k0/(1-c**2) + numpy.sqrt( (m*lambda0*k0/(1-c**2))**2 - \
        #                     ( (m * lambda0 * k0 / (1-c**2))**2 -1) )
        #
        # # Luca
        # sinalpha_luca = (-m * k0 * lambda0 / (c ** 2 - 1)) + \
        #             numpy.sqrt(1 + (m * m * c * c * k0 * k0 * lambda0 * lambda0) / (
        #                         (c ** 2 - 1) ** 2))
        #
        # print("sin_alpha numeric, me, luca: ",sinalpha,sinalpha_me,sinalpha_luca)

        # change sign of beta
        # sinbeta *= -1

        # alpha = numpy.arcsin(sinalpha)

        sinbeta = order * wavelength * line_density - sinalpha

        # beta = numpy.arcsin(sinbeta)

    return sinalpha,sinbeta

def vls_coefficients_calculate(sinalpha, sinbeta, rg, rgp, line_density=1000, wavelength=10e-10, order=1):

    #
    # calculate grating coefficients:
    #

    denominator = 2.0 * order * wavelength * line_density

    # VLS ShadowOui preprocessor
    # self.b2 = (((numpy.cos(alpha) ** 2) / self.r_a) + ((numpy.cos(beta) ** 2) / self.r_b)) / (
    #             -2 * m * self.k * wavelength)
    # self.b3 = ((numpy.sin(alpha) * numpy.cos(alpha) ** 2) / self.r_a ** 2 - \
    #            (numpy.sin(beta) * numpy.cos(beta) ** 2) / self.r_b ** 2) / (-2 * m * self.k * wavelength)
    # self.b4 = (((4 * numpy.sin(alpha) ** 2 - numpy.cos(alpha) ** 2) * numpy.cos(alpha) ** 2) / self.r_a ** 3 + \
    #            ((4 * numpy.sin(beta) ** 2 - numpy.cos(beta) ** 2) * numpy.cos(beta) ** 2) / self.r_b ** 3) / (
    #                       -8 * m * self.k * wavelength)


    b2 = ( (1-sinalpha**2)/rg + (1-sinbeta**2)/rgp) / denominator
    b3 = ( sinalpha*(1-sinalpha**2)/rg**2 + sinbeta*(1-sinbeta**2)/rgp**2) / denominator
    b4 = (((4 * sinalpha ** 2 - (1-sinalpha**2)) * (1-sinalpha**2)) / rg ** 3 + \
               ((4 * sinbeta ** 2 - (1-sinbeta**2)) * (1-sinbeta**2)) / rgp ** 3) / \
                (4 * denominator )

    # print(" b2=%f\n b3=%f\n b4=%f\n"%(b2,b3,b4))
    # print("Shadow coefficients:\n c1=%f\n c2=%f\n c3=%f\n"%(k0,2*b2*k0,-3*b3*k0))
    # shadow_coefficients = (k0, 2 * b2 * k0, -3 * b3 * k0, 4 * b4 * k0)
    # print("Shadow coefficients:\n c1=%f\n c2=%f\n c3=%f\n c4=%f\n" % shadow_coefficients)

    return b2,b3,b4

def vls_coefficients_convert_to_shadow(k0,b2,b3,b4):
    return k0, 2 * b2 * k0, -3 * b3 * k0, 4 * b4 * k0


if __name__ == "__main__":
    #
    # inputs
    #

    Sx = 23.7e-6
    Sy = 24.3e-6

    Sxp = 16.9e-6
    Syp = 16.7e-6


    L = 50.0

    m = 1 # order

    size_at_image_fwhm = 5e-6


    #
    # magnification and mirror and grating positions
    #

    Mm = size_at_image_fwhm / (2.35 * Sx)

    # mirror
    s = L*(1-Mm)
    sp = L*Mm

    sp = 4.12
    s = 45.88

    # grating
    rgp = sp + 1
    rg = L - rgp

    Mg = rgp / rg

    print("L=%f m, Mm=%f"%(L,Mm))
    print("M3 p=%f m, q=%f m, Mm=%f"%(s,sp,Mm))
    print("G p=%f m, q=%f m, Mg=%f"%(rg,rgp,Mg))

    #
    # solve for c (Eq. 5.6 ALS-U CDR)
    #
    c = rgp * s / (L * sp - rgp * sp)


    print("c = %f"%c)


    #
    # get angles
    #

    k0 = 114000.0
    energy0 = 603.0
    lambda0 = m2ev / energy0

    sinalpha,sinbeta = solve_grating_equation(line_density=k0,wavelength=lambda0,c=c,order=m)

    alpha = numpy.arcsin(sinalpha)
    beta = numpy.arcsin(sinbeta)
    theta = (alpha-beta)/2

    print("\nMain energy=%f eV, wavelength=%f A"%(energy0,lambda0*1e10))
    print("\nalpha=%f deg, beta=%f deg, theta=%f deg"%(
        alpha*180/numpy.pi,beta*180/numpy.pi,theta*180/numpy.pi))


    #
    # get resolving power
    #

    delta_lambda0_exit = k0**(-1) / numpy.abs(m) * numpy.cos(beta) * size_at_image_fwhm / rgp
    delta_lambda0_source = k0**(-1) / numpy.abs(m) * numpy.cos(alpha) * (2.35 * Sy) / rg
    delta_lambda0 = numpy.sqrt( delta_lambda0_source**2 + delta_lambda0_exit**2)

    print("Delta_lambda: %g, %g, %g"%(delta_lambda0_source,delta_lambda0_exit,delta_lambda0))
    print("Resolving power: %g, %g, %g"%(lambda0/delta_lambda0_source,lambda0/delta_lambda0_exit,lambda0/delta_lambda0))


    #
    # calculate grating coefficients:
    #

    b2,b3,b4 = vls_coefficients_calculate(sinalpha, sinbeta, rg, rgp, line_density=k0, wavelength=lambda0, order=m)

    print(" b2=%f\n b3=%f\n b4=%f\n"%(b2,b3,b4))


    sh0,sh1,sh2,sh3 = vls_coefficients_convert_to_shadow(k0,b2,b3,b4)

    print("VLS coefficients (SI):\n k0=%f\n b2=%f\n b3=%f\n b4=%f\n"%(k0,b2,b3,b4))
    print("VLS Shadow coefficients (SI):\n c1=%f\n c2=%f\n c3=%f\n c4=%f\n" % (sh0,sh1,sh2,sh3))