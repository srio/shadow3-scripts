
import numpy

from scipy import interpolate
from scipy import spatial

import matplotlib.pylab as plt
import matplotlib as mpl

# from scipy.spatial import Delaunay
# from scipy.interpolate import griddata


def load_comsol_file(filename,nx=300,ny=100,xmax=None,ymax=None,four_quadrants=2,do_plot=False):


    a = numpy.loadtxt(filename,skiprows=9)

    node  = numpy.round(a,10) # 1/4 model: 55mm*35mmm, units in m

    # Coordinates

    X0 = node[:,2]  # X=x in m
    Y0 = node[:,0]  # Y=z in m
    Z0 = node[:,3]  # Z=uy vertical displacement in m

    if four_quadrants == 2:
        X1 = numpy.concatenate( (X0,-X0) )
        Y1 = numpy.concatenate( (Y0, Y0) )
        Z1 = numpy.concatenate( (Z0, Z0) )

        X = numpy.concatenate( (X1, X1) )
        Y = numpy.concatenate( (Y1,-Y1) )
        Z = numpy.concatenate( (Z1, Z1) )
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
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()


    if four_quadrants == 2:
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


    if four_quadrants == 1:
        z1 = numpy.vstack((numpy.flip(z,axis=0),z[1:,:]))
        z2 = numpy.hstack((numpy.flip(z1,axis=1),z1[:,1:]))

        xlin2 = numpy.concatenate((-xlin[::-1],xlin[1:]))
        ylin2 = numpy.concatenate((-ylin[::-1],ylin[1:]))

        return z2,xlin2,ylin2
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



if __name__ == "__main__":

    from srxraylib.plot.gol import plot_image, plot

    z,x,y = load_comsol_file("brut.txt",nx=400,ny=150,four_quadrants=2,do_plot=False)

    plot_image(z,x,y,xtitle="X (%d pixels, max:%f)"%(x.size,x.max()),ytitle="Y (%d pixels, max:%f)"%(y.size,y.max()),show=0)
    plot(x,z[:,z.shape[1]//2],show=0)
    plot(y,z[z.shape[0]//2,:])

    print(z.shape,x.shape,y.shape)

    write_shadow_surface(z,x,y, filename="brut.dat")





