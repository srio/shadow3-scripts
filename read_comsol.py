
import numpy
from scipy.interpolate import interp2d


def load_comsol_file(file_name,points_in_x=50,points_in_y=250):

    if points_in_x %2 == 0:
        raise Exception("Please use an odd number of points in X")
    if points_in_y %2 == 0:
        raise Exception("Please use an odd number of points in Y")

    a  =numpy.loadtxt(file_name,skiprows=9)

    print("Mean   StDev  Min  Max")
    for i in range(4):
        print("%f  %g   %f   %f  "%(a[:,i].mean(),a[:,i].std(),a[:,i].min(),a[:,i].max()))

    x_c = a[:,0]
    y_c = a[:,1] + a[:,3]
    z_c = a[:,2]

    x_s = z_c
    y_s = x_c
    z_s = y_c


    f = interp2d(x_s, y_s, z_s, kind='linear')

    xnew = numpy.linspace(0.0,x_s.max(),int(1+points_in_x/2))
    ynew = numpy.linspace(0.0,y_s.max(),int(1+points_in_y/2))

    X = numpy.outer(xnew,numpy.zeros_like(ynew))
    Y = numpy.outer(numpy.zeros_like(xnew),ynew)

    Z = f(xnew,ynew).T


    XNEW = numpy.hstack((-xnew[::-1],xnew[1:]))
    YNEW = numpy.hstack((-ynew[::-1],ynew[1:]))

    ZNEW1 = numpy.zeros( (XNEW.size, YNEW.size))
    ZNEW = numpy.zeros( (XNEW.size, YNEW.size))

    for i in range(xnew.size):
        ZNEW1[i,:] = numpy.hstack((Z[i,::-1],Z[i,1:]))
    for j in range(YNEW.size):
        a1 = ZNEW1[0:xnew.size,j].copy()
        a1 = a1[::-1]
        a2 = ZNEW1[1:xnew.size,j]
        ZNEW[:,j] = numpy.hstack((a1,a2))

    return ZNEW,XNEW,YNEW

def write_shadow_surface(s,xx,yy,outFile='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,outFile=outFile)
      INPUTS:
           z - 2D array of heights
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           out - 1=Success, 0=Failure
           outFile - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """
    out = 1

    try:
       fs = open(outFile, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+outFile)
       return
    else:
        # dimensions
        fs.write( repr(xx.size)+" "+repr(yy.size)+" \n" )
        # y array
        for i in range(yy.size):
            fs.write(' ' + repr(yy[i]) )
        fs.write("\n")
        # for each x element, the x value and the corresponding z(y)
        # profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps = tmps + "  " + repr(s[i,j])
            fs.write(' ' + repr(xx[i]) + " " + tmps )
            fs.write("\n")
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+outFile+" written to disk.")



from srxraylib.plot.gol import plot, plot_image, plot_contour

Z,X,Y = load_comsol_file("/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/surf3_comsol.txt",points_in_x=51,points_in_y=551)
Z = Z - Z[int(X.size/2),int(Y.size/2)]

print(Z.shape,X.shape,Y.shape)
print(X)
plot_image(Z,X,Y,aspect='auto',show=False)
plot_contour(Z,X,Y)
write_shadow_surface(Z,X,Y)
print(X.max(),Y.max())
import os
print(os.getcwd())
plot(Y,1e6*Z[25,:],xtitle="Y [m]",ytitle="Z [um]")
