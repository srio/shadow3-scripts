
import numpy
from scipy import interpolate


# class GFile(object):
#
#     def __init__(self):
#
#         self.source = None
#         self.oe_list = []


if __name__ == "__main__":

    filename = "brut.txt"

    a = numpy.loadtxt(filename,skiprows=9)

    print(a.shape)

    node  = numpy.round(a,10) # 1/4 model: 55mm*35mmm, units in m

    # Changes coordinates

    X = 1e3 * node[:,0]  # X=x in mm
    Y = 1e3 * node[:,2]  # Y=z in mm
    Z = 1e6 * node[:,3]  # Z=uy vertical displacement in microns

    print(X.shape,Y.shape,Z.shape)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%% Interpolate: 3 methods
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #
    # %%%%% Grid for interpolation
    #
    xlin = numpy.linspace(0,55,300)
    ylin = numpy.linspace(0,35,100)
    #
    # [x,y]=meshgrid(xlin,ylin);

    x = numpy.outer(xlin,numpy.ones_like(ylin))
    y = numpy.outer(numpy.ones_like(xlin),ylin)

    print(x.shape,y.shape)

    f = interpolate.interp2d(X,Y,Z, kind='linear')

    z = f(xlin,ylin).T

    print(xlin.shape,ylin.shape,z.shape)

    from srxraylib.plot.gol import plot_image


    z1 = numpy.vstack((numpy.flip(z,axis=0),z[1:,:]))

    # plot_image(z1) #,xlin,ylin)

    print(z1.shape)

    z2 = numpy.hstack((numpy.flip(z1,axis=1),z1[:,1:]))



    print(z2.shape)

    x2 = numpy.concatenate((-xlin[::-1],xlin[1:]))
    print(x2)

    y2 = numpy.concatenate((-ylin[::-1],ylin[1:]))

    print(z2.shape,x2.shape,y2.shape)

    plot_image(z2,x2,y2)






