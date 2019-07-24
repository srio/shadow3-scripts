import Shadow
import numpy
from srxraylib.plot.gol import plot, plot_scatter
import matplotlib.pylab as plt

# ;+
# ;
# ;  NAME:
# ; 	RESPOWER
# ;  PURPOSE:
# ; 	to compute the resolving power E/DE from the dispersion in a plane
# ;       (usually the exit slit plane).
# ;  CATEGORY:
# ;        SHADOW's utilities
# ;  CALLING SEQUENCE:
# ; 	respower,'myfile',col_E, col_D
# ;  INPUTS:
# ; 	myfile  name of the file with data (between quotes)
# ;                (it can also be a structure like idl_var)
# ;       col_E:  the SHADOW column with the energy information:
# ;                   11: Energy in eV
# ;                   19: Wavelength in A
# ;       col_D:  the SHADOW column with the legth in the dispersion direction (1 or 3)
# ;  KEYWORD PARAMETERS:
# ; 		ERANGE=[emin,emax] range of the energy variable for the TOP plot
# ; 		ZRANGE=[Zmin,Zmax]              Z
# ; 		E2RANGE=[emin,emax] range of the energy variable for the BOTTOM plot
# ; 		Z2RANGE=[Zmin,Zmax]              Z
# ;               NBINS              number of bins for the histogram
# ; 		HLIMIT             the normalized heigh at which the histogram
# ;                                  limits taken for calculating the resolving power.
# ;                                  Default: 0.1
# ; 		TITLE='top_title' title to be written at the top
# ; 		NOPLOT = When set, inhibits the graph.
# ;
# ;  OUTPUTS:
# ; 	a plot
# ;  OPTIONAL OUTPUT PARAMETERS:
# ; 	resolvingPower: Set this keyword to a named variable to receive the Resolving Power
# ;       deltaEmin: Set this keyword to a named variable to receive the minimum DE
# ;                  (at zero exist slit opening)
# ;       Pendent: Set this keyword to a named variable to receive P
# ;
# ;                Note that for any exit slit aperture  DE = deltaE + DZ/P
# ;  COMMON BLOCKS:
# ; 	None.
# ;  SIDE EFFECTS:
# ; 	None.
# ;  RESTRICTIONS:
# ; 	None.
# ;  PROCEDURE:
# ; 	Compute the dispersion plot (Z vs DE, top graph) then calculate the dispersion
# ;       band by
# ;         i) Removing the regression line (bottom: plot with regression removed)
# ;        ii) Calculate Z histogram
# ;       iii) Calculate the histogram limits at HLIMIT height (bottom plot),
# ;        iv) Plot these limits in the upper graphic to define the dispersion band,
# ;         v) compute related parameters:
# ;          The "pendent" P from the regression fit gives the disperion DZ/DE
# ;          The "origin" is the Eo value at the upper boundary of the dispersion band
# ;          The Resolving Power is  Eo/P
# ;          The "Delta" value is the DE corresponding to the histogram limits, or
# ;              equivalently, the DE at zero slit opening
# ;
# ;          The energy bandwith at a given slit opening Z is DE= Delta+ Z/P
# ;
# ;  MODIFICATION HISTORY:
# ; 	by M. Sanchez del Rio. ESRF. Grenoble, around 1995
# ; 	2013/03/18 srio@esrf.eu documented, extracted some parameters
# ;
# ;-

# Pro respower,input,col1,col2,nolost=nolost,nbins=nbins, hlimit = hlimit, title=title, $
#    erange=erange, zrange=zrange,z2range=z2range, $
#    resolvingPower=resolvingPower,deltaEmin=deltaE,pendent=pendent

def respower(beam0,colE,col1,nolost=True,nbins=100,hlimit=0.1,do_plot=True,
             title="",erange=None, zrange=None,z2range=None,):


    if  colE != 11 and col1 != 19:
        raise Exception('RESPOWER: Warning: first column is NOT energy or wavelength.')

    beam = beam0.duplicate()

    #
    # get data
    #
    energy = beam.getshonecol(colE,nolost=nolost)
    z = beam.getshonecol(col1, nolost=nolost)


    degree=1
    coeff = numpy.polyfit(energy, z, degree, rcond=None, full=False, w=None, cov=False)

    yfit = coeff[1] + coeff[0] * energy
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
    print('The linear fit parameters are y = %f + %f '%(coeff[1],coeff[0]))
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


    if False:
        f = plot_scatter(beam.getshonecol(11), beam.getshonecol(col1),show=0)
        f[1].plot(energy, energy*0+coordinates_at_hlimit[0])
        f[1].plot(energy, energy*0+coordinates_at_hlimit[1])
        plt.show()

    if do_plot:
        if colE == 11:
            xtitle = "Photon energy [eV]"
        elif colE == 19:
            xtitle = "Photon wavelength [A]"

        ytitle = "column %i [user units]"%col1
        f = plot_scatter(beam.getshonecol(11), z,show=0,xtitle=xtitle,ytitle=ytitle,title="")
        f[1].plot(energy, yfit)
        f[1].plot(energy, yfit+coordinates_at_hlimit[0])
        f[1].plot(energy, yfit+coordinates_at_hlimit[1])
        f[1].plot(energy, energy*0)
        # IF (col1 EQ 19) THEN BEGIN
        #  oplot,[orig+deltax1,orig+deltax1],[-1000,1000]
        # ENDIF ELSE BEGIN
        #  oplot,[orig-deltax1,orig-deltax1],[-1000,1000]
        # ENDELSE
        if colE == 19: # wavelength
            f[1].plot(numpy.array((orig + deltax1, orig + deltax1)), numpy.array((-1000, 1000)))
            f[1].plot(numpy.array((orig - deltax2, orig - deltax2)), numpy.array((-1000, 1000)))
        else: # energy
            f[1].plot(numpy.array((orig - deltax1, orig - deltax1)), numpy.array((-1000, 1000)))
            f[1].plot(numpy.array((orig + deltax2, orig + deltax2)), numpy.array((-1000, 1000)))


        plt.show()

if __name__ == "__main__":

    beam0 = Shadow.Beam()
    beam0.load('/Users/srio/Oasys/star.04')
    respower(beam0,11,1)

