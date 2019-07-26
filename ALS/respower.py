import Shadow
import numpy


#+
#
#  NAME:
# 	respower
#  PURPOSE:
# 	to compute the resolving power E/DE from the dispersion in a plane (usually the exit slit plane).
#  CATEGORY:
#        SHADOW's utilities
#  CALLING SEQUENCE:
# 	respower(beam,col_E, col_D)
#  INPUTS:
# 	beam  a Shadow.Beam() object
#       col_E:  the SHADOW column with the energy information, either:
#                   11: Energy in eV
#                   19: Wavelength in A
#       col_D:  the SHADOW column with the legth in the dispersion direction (1 or 3)
#  KEYWORD PARAMETERS:
#       nbins              number of bins for the histogram
# 		hlimit             the normalized height at which the histogram limits are taken for calculating the
#                          resolving power. Default: 0.1
# 		title              title to be written at the top
# 		do_plot            True: display plots, False: No plots             Z
#
#  OUTPUTS:
#       a dictionaire with tags:
#       "resolvingPower"
#       "deltaE": the minimum DE (at zero exist slit opening)
#       "pendent":  P
#       "origin"
#
#       Note that for any exit slit aperture  DE = deltaE + DZ/P
#
#  PROCEDURE:
# 	Compute the dispersion plot (Z vs DE) then calculate the dispersion band by
#         i) Removing the regression line
#        ii) Calculate Z histogram
#       iii) Calculate the histogram limits at hlimit height (bottom plot),
#        iv) Plot these limits in the dispersion graphic (Z vs E) to define the dispersion band,
#         v) compute related parameters:
#          The "pendent" P from the regression fit gives the disperion DZ/DE
#          The "origin" is the Eo value at the upper boundary of the dispersion band
#          The Resolving Power is  Eo/P
#          The "Delta" value is the DE corresponding to the histogram limits, or
#              equivalently, the DE at zero slit opening
#
#          The energy bandwith at a given slit opening Z is DE= Delta+ Z/P
#
#  MODIFICATION HISTORY:
# 	by M. Sanchez del Rio. ESRF. Grenoble, around 1995
# 	2013/03/18 srio@esrf.eu documented, extracted some parameters
#   2019/07/24 srio@lbl.gov python version
#
#-

def respower(beam0,colE,col1,nolost=True,nbins=100,hlimit=0.1,do_plot=True,title=""):


    if  colE == 11 or colE == 19:
        pass
    else:
        raise Exception('First column is NOT energy or wavelength.')

    beam = beam0.duplicate()

    #
    # get data
    #
    energy = beam.getshonecol(colE,nolost=nolost)
    energy_all_rays = beam.getshonecol(colE, nolost=False)
    z = beam.getshonecol(col1, nolost=nolost)


    degree=1
    coeff = numpy.polyfit(energy, z, degree, rcond=None, full=False, w=None, cov=False)

    yfit = coeff[1] + coeff[0] * energy_all_rays
    beam.rays[:,col1-1] -= yfit

    #
    # histogram
    #
    tkt = beam.histo1(col1,nbins=nbins,nolost=nolost,calculate_widths=True)
    hx = tkt["bin_center"]
    hy = tkt["histogram"]

    # ;
    # ; get histo maximum and limits
    # ;
    tt = numpy.where(hy >= (hy.max() * hlimit))
    if hy[tt].size > 1:
        binSize = hx[1] - hx[0]
        width_at_hlimit = binSize * (tt[0][-1] - tt[0][0])
        coordinates_at_hlimit = (hx[tt[0][0]], hx[tt[0][-1]])
        coordinates_at_center = hx[numpy.argmax(hy)]

    else:
        raise Exception("Failed to compute width at %f height"%hlimit)

    deltax1 = numpy.abs((coordinates_at_hlimit[0] - coordinates_at_center) / coeff[0])
    deltax2 = numpy.abs((coordinates_at_hlimit[1] - coordinates_at_center) / coeff[0])
    deltaE = deltax1 + deltax2
    # deltaE = numpy.abs((coordinates_at_hlimit[1] - coordinates_at_hlimit[0]) / coeff[0])


    # ;
    # ; data
    # ;
    orig = -1.0 * (coeff[1] + coordinates_at_center) / coeff[0]
    resolvingPower = orig/deltaE
    pendent = coeff[0]

    print("\n\n*********************************************************")
    print('The linear fit parameters are y = %f + %f x'%(coeff[1],coeff[0]))
    print('Mean of residuals: %g'%beam.rays[:,col1-1].mean())
    print('StDev of residuals: %g'%beam.rays[:,col1-1].std())
    print("width_at_hlimit: %g "%width_at_hlimit)
    print("coordinates_at_hlimit = (%g,%g)" % (coordinates_at_hlimit))
    print('Linear fit pendent P: %f'%pendent)
    print('Histogram peak at: %g'%coordinates_at_center)
    print('Histogram base line: %f'%hlimit)
    print('DeltaE(DZ=0) = %f'%deltaE)
    print('Origin E = %f'%orig)
    print('Resolving Power = %f'%resolvingPower)
    print('Intensity = %f'%beam.intensity(nolost=nolost))
    print('DE ~ DeltaE + DZ/|P|')
    print("*********************************************************\n\n")

    return {"colE":colE,"col1":col1,"nbins":nbins,"nolost":nolost,"hlimit":hlimit,"title":title, # inputs
            "resolvingPower":resolvingPower,
            "origin":orig,
            "pendent":pendent,
            "deltaE":deltaE,
            "coeff":coeff,
            "coordinates_at_hlimit":coordinates_at_hlimit,
            "coordinates_at_center":coordinates_at_center,
            "deltax1":deltax1,
            "deltax2":deltax2,
            "histo_dict":tkt}

def respower_plot(beam,d,plot_substracted=False,nolost=True):
    from srxraylib.plot.gol import plot, plot_scatter
    import matplotlib.pylab as plt

    colE = d["colE"]
    col1 = d["col1"]
    coeff = d["coeff"]
    nolost = d["nolost"]
    coordinates_at_hlimit = d["coordinates_at_hlimit"]
    orig = d["origin"]
    title = d["title"]
    deltax1 = d["deltax1"]
    deltax2 = d["deltax2"]


    if colE == 11:
        xtitle = "Photon energy [eV]"
        unit = "eV"
    elif colE == 19:
        xtitle = "Photon wavelength [A]"
        unit = "A"

    ytitle = "column %i [user units]"%col1

    energy = beam.getshonecol(colE,nolost=nolost)
    z = beam.getshonecol(col1, nolost=nolost)
    yfit = coeff[1] + coeff[0] * energy

    #
    # substracted plot
    #
    if plot_substracted:
        f = plot_scatter(energy, z-(coeff[1]+coeff[0]*energy),xtitle=xtitle, ytitle=ytitle, title=title,show=0)
        f[1].plot(energy, energy*0+coordinates_at_hlimit[0])
        f[1].plot(energy, energy*0+coordinates_at_hlimit[1])
        plt.show()

    #
    # main plot
    #

    g = plot_scatter(energy, z,show=0,xtitle=xtitle,ytitle=ytitle,
                     title=title+" E/DE=%d, DE=%f %s"%(d["resolvingPower"],d["deltaE"],unit))
    g[1].plot(energy, yfit)
    g[1].plot(energy, yfit+coordinates_at_hlimit[0])
    g[1].plot(energy, yfit+coordinates_at_hlimit[1])
    g[1].plot(energy, energy*0)
    if colE == 19: # wavelength
        g[1].plot(numpy.array((orig + deltax1, orig + deltax1)), numpy.array((-1000, 1000)))
        g[1].plot(numpy.array((orig - deltax2, orig - deltax2)), numpy.array((-1000, 1000)))
    else: # energy
        g[1].plot(numpy.array((orig - deltax1, orig - deltax1)), numpy.array((-1000, 1000)))
        g[1].plot(numpy.array((orig + deltax2, orig + deltax2)), numpy.array((-1000, 1000)))


    plt.show()


if __name__ == "__main__":

    beam0 = Shadow.Beam()
    # file = '/Users/srio/Oasys/star.04'

    file = "C:\\Users\\Manuel\\Oasys\\star_slit.dat"

    beam0.load(file)



    dict = respower(beam0,11,1,hlimit=0.5,nolost=True)
    for key in dict.keys():
        print(key," = ",dict[key])

    respower_plot(beam0,dict,nolost=True)



