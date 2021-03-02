#
# script to make the calculations (created by XOPPY:undulator_spectrum)
#
def run_power_density():
    from orangecontrib.xoppy.util.xoppy_undulators import xoppy_calc_undulator_power_density

    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"] = 6.0
    h5_parameters["ELECTRONENERGYSPREAD"] = 0.00138
    h5_parameters["ELECTRONCURRENT"] = 0.2
    h5_parameters["ELECTRONBEAMSIZEH"] = 1.48e-05         * 10
    h5_parameters["ELECTRONBEAMSIZEV"] = 3.7e-06          * 10
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 2.8e-06    * 10
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 1.5e-06    * 10
    h5_parameters["PERIODID"] = 0.025
    h5_parameters["NPERIODS"] = 184.0
    h5_parameters["KV"] = 1.862765
    h5_parameters["KH"] = 0.0
    h5_parameters["KPHASE"] = 0.0
    h5_parameters["DISTANCE"] = 28.0
    h5_parameters["GAPH"] = 0.01
    h5_parameters["GAPV"] = 0.01
    h5_parameters["HSLITPOINTS"] = 41 # 81
    h5_parameters["VSLITPOINTS"] = 41 # 81
    h5_parameters["METHOD"] = 2 # 1 is URGENT
    h5_parameters["USEEMITTANCES"] = 1
    h5_parameters["MASK_FLAG"] = 0
    h5_parameters["MASK_ROT_H_DEG"] = 0.0
    h5_parameters["MASK_ROT_V_DEG"] = 88.0
    h5_parameters["MASK_H_MIN"] = -50.0
    h5_parameters["MASK_H_MAX"] = 50.0
    h5_parameters["MASK_V_MIN"] = -5.0
    h5_parameters["MASK_V_MAX"] = 5.0

    h, v, p, code = xoppy_calc_undulator_power_density(
        ELECTRONENERGY=h5_parameters["ELECTRONENERGY"],
        ELECTRONENERGYSPREAD=h5_parameters["ELECTRONENERGYSPREAD"],
        ELECTRONCURRENT=h5_parameters["ELECTRONCURRENT"],
        ELECTRONBEAMSIZEH=h5_parameters["ELECTRONBEAMSIZEH"],
        ELECTRONBEAMSIZEV=h5_parameters["ELECTRONBEAMSIZEV"],
        ELECTRONBEAMDIVERGENCEH=h5_parameters["ELECTRONBEAMDIVERGENCEH"],
        ELECTRONBEAMDIVERGENCEV=h5_parameters["ELECTRONBEAMDIVERGENCEV"],
        PERIODID=h5_parameters["PERIODID"],
        NPERIODS=h5_parameters["NPERIODS"],
        KV=h5_parameters["KV"],
        KH=h5_parameters["KH"],
        KPHASE=h5_parameters["KPHASE"],
        DISTANCE=h5_parameters["DISTANCE"],
        GAPH=h5_parameters["GAPH"],
        GAPV=h5_parameters["GAPV"],
        HSLITPOINTS=h5_parameters["HSLITPOINTS"],
        VSLITPOINTS=h5_parameters["VSLITPOINTS"],
        METHOD=h5_parameters["METHOD"],
        USEEMITTANCES=h5_parameters["USEEMITTANCES"],
        MASK_FLAG=h5_parameters["MASK_FLAG"],
        MASK_ROT_H_DEG=h5_parameters["MASK_ROT_H_DEG"],
        MASK_ROT_V_DEG=h5_parameters["MASK_ROT_V_DEG"],
        MASK_H_MIN=h5_parameters["MASK_H_MIN"],
        MASK_H_MAX=h5_parameters["MASK_H_MAX"],
        MASK_V_MIN=h5_parameters["MASK_V_MIN"],
        MASK_V_MAX=h5_parameters["MASK_V_MAX"],
        h5_file="undulator_power_density.h5",
        h5_entry_name="XOPPY_POWERDENSITY",
        h5_initialize=True,
        h5_parameters=h5_parameters,
    )
    return p, h, v


if __name__ == "__main__":
    # example plot
    from srxraylib.plot.gol import plot_image
    import numpy

    p, h, v = run_power_density()
    plot_image(p, h, v, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2")


    from syned.storage_ring.magnetic_structures.undulator import Undulator
    from syned.storage_ring.electron_beam import ElectronBeam

    ebeam = ElectronBeam(energy_in_GeV=6.0, current=0.2)
    K_vertical = 1.862765
    und = Undulator(K_vertical=K_vertical, period_length=0.025, number_of_periods=184)

    gamma = ebeam.gamma()
    distance = 28.0

    H = numpy.outer(h, numpy.ones_like(v))
    V = numpy.outer(numpy.ones_like(h), v)
    theta = numpy.sqrt(H**2 + V**2) / (distance * 1000)
    phi = numpy.arctan2(V, H)

    alpha = gamma * theta

    n = 1.0 # harmonic
    xi = n / (1 + 0.5 * K_vertical**2 + alpha**2)

    X = 2 * xi * alpha * (K_vertical * numpy.cos(phi))
    Y = xi * K_vertical**4 / 4
    Phi = 0.0

    # plot_image(X, title="X")
    # plot_image(Y, title="Y")

    import scipy.special as special

    S1 = numpy.zeros_like(X)
    S2 = numpy.zeros_like(X)

    for i in range(-30,30,1):
        JNX = special.jn(i,X)
        JNY = special.jn(i,Y)
        S1 += special.jn(i,Y) * special.jn(2*i+n,X)
        S1 += i * special.jn(i,Y) * special.jn(2*i+n,X)
        print(">>>>", i)

    print(numpy.abs(JNX).max(), numpy.abs(JNY).max(), )

    Ax = xi * (2 * alpha * numpy.cos(phi) * S1 - K_vertical * (2 * n * S1 + 4 * S2) / X)
    Ay = xi * (2 * alpha * numpy.sin(phi) * S1 )

    A2 = Ax**2 + Ay**2
    PWR = A2 / (1 + K_vertical**2 + alpha**2)
    plot_image(A2, title="A")
    plot_image(PWR, title="PWR")
    # plot_image(S2, title="S2")

