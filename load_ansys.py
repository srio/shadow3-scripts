import numpy

from scipy import interpolate
from scipy import spatial

import matplotlib.pylab as plt
import matplotlib as mpl


def load_ansys_file(filename, nx=300, ny=100, xmax=None, ymax=None, mirror_data=2, do_plot=False,to_meter=1e-3):


    a = numpy.loadtxt(filename,skiprows=9)

    node  = numpy.round(a,10)

    # Coordinates

    X0 = node[:,1] * to_meter  # X=x in m
    Y0 = node[:,2] * to_meter  # Y=z in m
    Z0 = node[:,4] * to_meter  # Z=uy vertical displacement in m

    if mirror_data == 2: # mirror on X and Y (quadrant)
        X1 = numpy.concatenate( (X0,-X0) )
        Y1 = numpy.concatenate( (Y0, Y0) )
        Z1 = numpy.concatenate( (Z0, Z0) )

        X = numpy.concatenate( (X1, X1) )
        Y = numpy.concatenate( (Y1,-Y1) )
        Z = numpy.concatenate( (Z1, Z1) )
    elif mirror_data == 1:  # mirror on Y
        X = numpy.concatenate( (X0, X0) )
        Y = numpy.concatenate( (Y0,-Y0) )
        Z = numpy.concatenate( (Z0, Z0) )
    else:
        X = X0
        Y = Y0
        Z = Z0

    #
    # Interpolation
    #
    if xmax is None:
        xmax = numpy.abs(X).max()
    if ymax is None:
        ymax = numpy.abs(Y).max()


    # triangulation

    Pi = numpy.array([X, Y]).transpose()
    tri = spatial.Delaunay(Pi)

    if do_plot:
        plt.triplot(X0, Y0 , tri.simplices.copy())
        plt.plot(X0, Y0, "or", label = "Data")
        plt.grid()
        plt.legend()
        plt.xlabel("x TRI")
        plt.ylabel("y TRI")
        plt.show()


    if mirror_data == 2:
        xlin = numpy.linspace(-xmax,xmax,nx)
        ylin = numpy.linspace(-ymax,ymax,ny)
    elif mirror_data == 1:
        xlin = numpy.linspace(-xmax,xmax,nx)
        ylin = numpy.linspace(-ymax,ymax,ny)
    else:
        xlin = numpy.linspace(0,xmax,nx)
        ylin = numpy.linspace(0,ymax,ny)

    XLIN =  numpy.outer(xlin,numpy.ones_like(ylin))
    YLIN =  numpy.outer(numpy.ones_like(xlin),ylin)


    # interpolation
    P = numpy.array([XLIN.flatten(), YLIN.flatten() ]).transpose()

    z = interpolate.griddata(Pi, Z, P, method = "cubic").reshape([nx,ny])

    if do_plot:
        plt.contourf(XLIN, YLIN, z, 50, cmap = mpl.cm.jet)
        plt.colorbar()
        plt.contour(XLIN, YLIN, z, 20, colors = "k")
        #plt.triplot(Xi, Yi , tri.simplices.copy(), color = "k")
        plt.plot(X0, Y0, "or", label = "Data")
        plt.legend()
        plt.grid()
        plt.show()


    if mirror_data == 2:
        z1 = numpy.vstack((numpy.flip(z,axis=0),z[1:,:]))
        z2 = numpy.hstack((numpy.flip(z1,axis=1),z1[:,1:]))

        xlin2 = numpy.concatenate((-xlin[::-1],xlin[1:]))
        ylin2 = numpy.concatenate((-ylin[::-1],ylin[1:]))

        return z2,xlin2,ylin2
    elif mirror_data == 1:
        return z, xlin, ylin
    else:
        return z,xlin,ylin

def write_shadow_surface(s,xx,yy,filename='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,filename=filename)
      INPUTS:
           z - 2D array of heights z(x,y)
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           filename - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """

    try:
       fs = open(filename, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+filename)
       return
    else:
        # dimensions
        fs.write( "%d  %d \n"%(xx.size,yy.size))
        # y array
        for i in range(yy.size):
            fs.write("%g  "%(yy[i]))
        fs.write("\n")
        # for each x element, the x value followed by the corresponding z(y) profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps += "%g  "%(s[i,j])
            fs.write("%g    %s \n"%(xx[i],tmps))
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+filename+" written to disk.")



from srxraylib.plot.gol import plot_image

filename = "/scisoft/data/srio/EBS_READINESS/ID27/surfDeform_ID27-ML-DMM-L300Wcut6.6_U18-K2.1-OASYSCS.txt"
z,x,y = load_ansys_file(filename, nx=200, ny=50, mirror_data=1, do_plot=False)

plot_image(z,x,y,xtitle="X (%d pixels, max:%f)"%(x.size,x.max()),ytitle="Y (%d pixels, max:%f)"%(y.size,y.max()),aspect="auto")

write_shadow_surface(z,x,y, filename="id27_ml1.dat")

out_object = "id27_ml1.dat"

