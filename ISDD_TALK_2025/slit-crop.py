import numpy
import matplotlib.pylab as plt

def crop_center(img,cropx,cropy):
    x,y = img.shape
    startx = int(x*(1/2-cropx/2))
    starty = int(y*(1/2-cropy/2))
    return img[startx:startx+int(x*cropx),starty:starty+int(y*cropy)]

def plot_image_contour(p, h, v, aspect='equal', title="", xtitle="", ytitle="", show=0,
                       xrange=None,yrange=None, contour_step=100.0):
    # f = plt.imshow(p,h,v) #,aspect=aspect,title=title,xtitle=xtitle,ytitle=ytitle)

    fig = plt.figure()

    # cmap = plt.cm.Greys
    plt.imshow(p.T,origin='lower',extent=[h[0],h[-1],v[0],v[-1]],cmap=None,aspect=aspect)
    if True:
        plt.colorbar()
    ax = fig.gca()
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

    plt.title(title)

    if xrange is not None:
        plt.xlim( xrange )
    if yrange is not None:
        plt.ylim( yrange )

    levels = numpy.arange(0.0, p.max()*0.95, contour_step)  # Boost the upper limit to avoid truncation errors.

    vv = plt.axis()
    ff = plt.contour(p.T, levels, colors='k', origin='lower', extent=[h[0],h[-1],v[0],v[-1]])
    plt.clabel(ff, fmt='%d', colors='b', fontsize=14)
    plt.axis(vv)

    if show:
        plt.show()


def calculate_wiggler():
    #
    # script to make the calculations (created by XOPPY:undulator_spectrum)
    #
    from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_power_density

    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"] = 6.0
    h5_parameters["ELECTRONENERGYSPREAD"] = 0.001
    h5_parameters["ELECTRONCURRENT"] = 0.2
    h5_parameters["ELECTRONBEAMSIZEH"] = 3.12305e-05
    h5_parameters["ELECTRONBEAMSIZEV"] = 5.15684e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 4.51388e-06
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 1.93917e-06
    h5_parameters["PERIODID"] = 0.15
    h5_parameters["NPERIODS"] = 10.667
    h5_parameters["KV"] = 22.591
    h5_parameters["KH"] = 0.0
    h5_parameters["KPHASE"] = 0.0
    h5_parameters["DISTANCE"] = 23.7
    h5_parameters["GAPH"] = 0.1
    h5_parameters["GAPV"] = 0.1
    h5_parameters["HSLITPOINTS"] = 341
    h5_parameters["VSLITPOINTS"] = 341
    h5_parameters["METHOD"] = 2
    h5_parameters["USEEMITTANCES"] = 1
    h5_parameters["MASK_FLAG"] = 0
    h5_parameters["MASK_ROT_H_DEG"] = 0.0
    h5_parameters["MASK_ROT_V_DEG"] = 0.0
    h5_parameters["MASK_H_MIN"] = -1000.0
    h5_parameters["MASK_H_MAX"] = 1000.0
    h5_parameters["MASK_V_MIN"] = -1000.0
    h5_parameters["MASK_V_MAX"] = 1000.0

    horizontal, vertical, power_density, code = xoppy_calc_undulator_power_density(
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
    # example plot
    # from srxraylib.plot.gol import plot_image
    # plot_image(power_density, horizontal, vertical, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2")
    #
    # end script
    #
    return power_density, horizontal, vertical, h5_parameters["GAPH"], h5_parameters["DISTANCE"]


def calculate():
    #
    # script to make the calculations (created by XOPPY:undulator_spectrum)
    #
    from xoppylib.sources.xoppy_undulators import xoppy_calc_undulator_power_density

    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"] = 6.0
    h5_parameters["ELECTRONENERGYSPREAD"] = 0.001
    h5_parameters["ELECTRONCURRENT"] = 0.2
    h5_parameters["ELECTRONBEAMSIZEH"] = 0.0
    h5_parameters["ELECTRONBEAMSIZEV"] = 0.0
    h5_parameters["ELECTRONBEAMDIVERGENCEH"] = 0.0
    h5_parameters["ELECTRONBEAMDIVERGENCEV"] = 0.0
    h5_parameters["PERIODID"] = 0.018
    h5_parameters["NPERIODS"] = 111.111
    h5_parameters["KV"] = 1.4789
    h5_parameters["KH"] = 0.0
    h5_parameters["KPHASE"] = 0.0
    h5_parameters["DISTANCE"] = 23.7
    h5_parameters["GAPH"] = 0.01
    h5_parameters["GAPV"] = 0.01
    h5_parameters["HSLITPOINTS"] = 81
    h5_parameters["VSLITPOINTS"] = 81
    h5_parameters["METHOD"] = 2
    h5_parameters["USEEMITTANCES"] = 1
    h5_parameters["MASK_FLAG"] = 0
    h5_parameters["MASK_ROT_H_DEG"] = 0.0
    h5_parameters["MASK_ROT_V_DEG"] = 0.0
    h5_parameters["MASK_H_MIN"] = -1000.0
    h5_parameters["MASK_H_MAX"] = 1000.0
    h5_parameters["MASK_V_MIN"] = -1000.0
    h5_parameters["MASK_V_MAX"] = 1000.0

    horizontal, vertical, power_density, code = xoppy_calc_undulator_power_density(
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
    # example plot
    # from srxraylib.plot.gol import plot_image
    # plot_image(power_density, horizontal, vertical, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2")
    #
    # end script
    #

    return power_density, horizontal, vertical, h5_parameters["GAPH"], h5_parameters["DISTANCE"]

if __name__ == "__main__":

    #
    # PowerThroughSlit
    #
    CROP = numpy.linspace(1, 25, 50) / 25  # numpy.array([1.0,10/25,5/25,1/25])
    slit_aperture = numpy.zeros_like(CROP)
    slit_power = numpy.zeros_like(slit_aperture)

    # p, h, v, GAPH, DISTANCE = calculate()
    p, h, v, GAPH, DISTANCE = calculate_wiggler()



    for i, crop in enumerate(CROP):
        slit_aperture[i] = 1e3 * GAPH * crop
        slit_power[i] = crop_center(p, crop, crop).sum() * (h[1] - h[0]) * (v[1] - v[0])
        txt = "\n  Total power at at %5.2f m in %4.1f x %4.1f mm slit : %8.3f W" % (DISTANCE,
                                                                                    slit_aperture[i], slit_aperture[i],
                                                                                    slit_power[i])
        print(txt)

    #
    # PowerThroughSlit2D
    #
    hplus = numpy.array(numpy.where(h >= 0))
    hplus.shape = -1
    vplus = numpy.array(numpy.where(v >= 0))
    vplus.shape = -1
    p2 = p[hplus[0]:1 + hplus[-1], vplus[0]:1 + vplus[-1]].copy()
    p3 = numpy.cumsum(p2, 0)
    p3 = numpy.cumsum(p3, 1)
    p3 *= numpy.abs((h[1] - h[0]) * (v[1] - v[0])) * 4  # the four is because we double the interval in H and V

    hh = 2.0 * h[hplus]  # axes are now gap (aperture)
    vv = 2.0 * v[vplus]  # axes are now gap (aperture)


    from srxraylib.plot.gol import plot_image

    plot_image(p, h, v, xtitle="H [mm]", ytitle="V [mm]", title="Power density W/mm2", show=0)
    # plot_image(p3, hh, vv, title="Power [W] Through Slit", xtitle="X gap [mm]", ytitle="Y gap [mm]", show=0)
    plot_image_contour(p3, hh, vv, title="Power [W] Through Slit", xtitle="X gap [mm]", ytitle="Y gap [mm]", contour_step=1000, show=1)