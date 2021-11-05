import numpy

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
    #    i=0          i=1          i=2            i=3          i=4           i=5           i=6
    #j=0 1.178405e-08 1.156881e-13 -2.328567e-07 -3.826699e-09 -5.402403e+00 3.954034e-05 3.124443e+03
    #j=1 -8.743569e-08 -3.203787e-11 4.005674e-06 1.128360e-06 1.281728e+02 -8.909678e-03 -5.091596e+04
    #j=2 4.122065e-11 -1.747475e-08 9.937690e-01 4.738588e-04 -9.335430e+02 -3.952495e+00 3.619761e+05
    #j=3 3.015938e-08 1.103960e-06 -4.122661e+00 -3.998111e-02 3.596963e+03 3.424929e+02 -1.088331e+06
    #j=4 5.812141e-07 4.678890e-04 7.997619e+00 -1.170459e+01 -9.336506e+03 8.438191e+04 2.295149e+07
    #j=5 -2.233909e-04 -7.150309e-03 7.392859e+00 2.604628e+02 -2.897272e+04 -2.325954e+06 -2.010028e+08
    #j=6 -2.432996e-03 -2.834544e+00 -1.390798e+03 6.727502e+04 3.285300e+07 -4.479793e+08 -1.911958e+11


    coeffs = [
    [ 1.178405e-08 , 1.156881e-13 ,  -2.328567e-07, -3.826699e-09,  -5.402403e+00,  3.954034e-05,  3.124443e+03 ],
    [ -8.743569e-08, -3.203787e-11,   4.005674e-06,  1.128360e-06,   1.281728e+02, -8.909678e-03, -5.091596e+04 ],
    [ 4.122065e-11 , -1.747475e-08,   9.937690e-01,  4.738588e-04,  -9.335430e+02, -3.952495e+00,  3.619761e+05 ],
    [ 3.015938e-08 ,  1.103960e-06,  -4.122661e+00, -3.998111e-02,   3.596963e+03,  3.424929e+02, -1.088331e+06 ],
    [ 5.812141e-07 ,  4.678890e-04,   7.997619e+00, -1.170459e+01,  -9.336506e+03,  8.438191e+04,  2.295149e+07 ],
    [ -2.233909e-04, -7.150309e-03,   7.392859e+00,  2.604628e+02,  -2.897272e+04, -2.325954e+06, -2.010028e+08 ],
    [ -2.432996e-03, -2.834544e+00,  -1.390798e+03,  6.727502e+04,   3.285300e+07, -4.479793e+08, -1.911958e+11 ],
    ]

    coeffs = numpy.array(coeffs).T

    print("coeff 0,2: ", coeffs[0,2])

    nx = 201
    ny = 401
    x = numpy.linspace(-0.0125,0.0125,nx)
    y = numpy.linspace(-0.0125,0.0125,ny)


    X = numpy.outer(x, numpy.ones_like(y))
    Y = numpy.outer(numpy.ones_like(x), y)
    Z = numpy.zeros((x.size, y.size))

    for i in range(coeffs.shape[0]):
        for j in range(coeffs.shape[1]):
            Z += coeffs[i,j] * X**i * (-Y**j) # changed with respect Ken's eq.


    Zmean = (Z.max() - Z.min()) / 2
    print(Z.max() , Z.min(), Zmean)
    Z -= Z[nx//2,ny//2]
    # Z -= Z.min() + Zmean


    from srxraylib.plot.gol import plot_image

    plot_image(Z, x, y)


    write_shadow_surface(Z,x,y,filename='/Users/srio/Oasys/crystal_mesh.dat')



