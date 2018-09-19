
import numpy
from scipy import interpolate


def load_comsol_file(filename,nx=300,ny=100,xmax=None,ymax=None,kind='linear',four_quadrants=True):


    a = numpy.loadtxt(filename,skiprows=9)

    node  = numpy.round(a,10) # 1/4 model: 55mm*35mmm, units in m

    # Coordinates

    X = node[:,0]  # X=x in m
    Y = node[:,2]  # Y=z in m
    Z = node[:,3]  # Z=uy vertical displacement in m

    #
    # Interpolation
    #
    if xmax is None:
        xmax = X.max()
    if ymax is None:
        ymax = Y.max()

    xlin = numpy.linspace(0.0,xmax,nx)
    ylin = numpy.linspace(0.0,ymax,ny)


    f = interpolate.interp2d(X,Y,Z, kind=kind)

    z = f(xlin,ylin).T


    if four_quadrants:

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

    from srxraylib.plot.gol import plot_image

    z,x,y = load_comsol_file("brut.txt",nx=40,ny=150,four_quadrants=True,kind='linear')

    # plot_image(z,x,y,xtitle="X (%d pixels, max:%f)"%(x.size,x.max()),ytitle="Y (%d pixels, max:%f)"%(y.size,y.max()),)

    print(z.shape,x.shape,y.shape)

    write_shadow_surface(z,x,y, filename="brut.dat")





