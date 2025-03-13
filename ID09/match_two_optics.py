import numpy
from srxraylib.plot.gol import plot

def get_f2(ff1):
    q1 = 1 / (1 / ff1 - 1 / p1)
    p2 = D - q1
    ff2 = 1.0 / (1 / p2 + 1 / q2)
    return ff2, (q1 / p1) * (q2 / p2)

if __name__ == "__main__":
    p1 = 44.54
    q2 = 0.5
    D = 11.75 - q2



    F1 = numpy.linspace(5, 15, 100)
    F2, MM = get_f2(F1)

    M = (1 - q2 / F2) / (1 - p1 / F1)

    f2_asymp = 1/( 1/(p1 + D) + 1/(q2) )



    import xraylib
    radius = 200e-6
    symbol = "Be"
    density = 1.845
    photon_energy_ev = 15000.0
    delta = 1.0 - xraylib.Refractive_Index_Re(symbol, photon_energy_ev * 1e-3, density)
    print("delta: %g" % delta)
    print("Re(n): %g" % (1-delta))
    print("f_lens: %g" % (0.5 * radius / delta))

    R_torus = 2 * F1 / numpy.sin(2.49e-3)
    r_torus = 2 * F1 * numpy.sin(2.49e-3)
    for i in range(F1.size):
        if F2[i] > 0: print("f1: %5.3f, f2: %5.3f, M: %5.3f, N=%d, R_torus=%5.3f, r_torus=%5.3f" % (F1[i], F2[i], 1/numpy.abs(M[i]), radius/2/delta/F2[i], R_torus[i], r_torus[i]))


    #
    # plots
    #
    plot(F1,F2,
         F1, F2*0 + f2_asymp,
         xtitle='f1 [m]', ytitle='f2 [m]', title='(F1,F2)', yrange=[0,3], show=0)


    print("f2 asymp: ", f2_asymp)
    plot(
        F1, numpy.abs(M),
        F1, M*0 + (11.75/p1),
        xtitle='f1 [m]', ytitle='M',
        title="|M|", show=0)

    plot(
        F1, numpy.abs(1/M),
        F1, M*0 + (p1/11.75),
        xtitle='f1 [m]', ytitle='1/M',
        title="1/|M|")

    import xraylib
    symbol = "Be"
    density = 1.845
    photon_energy_ev = 15000.0
    delta = 1.0 - xraylib.Refractive_Index_Re(symbol, photon_energy_ev * 1e-3, density)
    print("delta: %g" % delta)