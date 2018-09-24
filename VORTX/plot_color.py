import numpy
import numpy as np
import h5py
from srxraylib.plot.gol import plot_image,plot_surface

# see https://stackoverflow.com/questions/17044052/mathplotlib-imshow-complex-2d-array

from colorsys import hls_to_rgb
from matplotlib.colors import hsv_to_rgb

import pylab as plt
from matplotlib.colors import Normalize



# def colorize(z):
#     n,m = z.shape
#     c = np.zeros((n,m,3))
#     c[np.isinf(z)] = (1.0, 1.0, 1.0)
#     c[np.isnan(z)] = (0.5, 0.5, 0.5)
#
#     idx = ~(np.isinf(z) + np.isnan(z))
#     A = np.angle(z[idx]) # (np.angle(z[idx]) + np.pi) / (2*np.pi)
#     # A = (A + 0.5) % 1.0
#     B = 1.0 - 1.0/(1.0+abs(z[idx])**0.3)
#     c[idx] = [hls_to_rgb(a, b, 0.8) for a,b in zip(A,B)]
#     return c
def Complex2HSV(z, rmin=None, rmax=None, hue_start=0):
    # get amplidude of z and limit to [rmin, rmax]
    amp = np.log10(np.abs(z)**2)
    amp -= amp.min()
    print("amp interval", amp.min(),amp.max())
    if rmin is not None:
        amp = np.where(amp < rmin, rmin, amp)
    if rmax is not None:
        amp = np.where(amp > rmax, rmax, amp)

    # ph = np.angle(z, deg=1) + hue_start
    # HSV are values in range [0,1]
    h = (np.angle(z, deg=1) %360) / 360 # (ph % 360) / 360
    s = 0.85 * np.ones_like(h)
    v = amp / amp.max()# (amp -rmin) / (rmax - rmin)
    return hsv_to_rgb(np.dstack((h,s,v)))

# def colorize(z):
#     r = np.abs(z)
#     arg = np.angle(z)
#
#     h = (arg + np.pi)  / (2 * np.pi) + 0.5
#     l = 1.0 - 1.0/(1.0 + r**0.3)
#     s = 0.8
#
#     c = np.vectorize(hls_to_rgb) (h,l,s) # --> tuple
#     c = np.array(c)  # -->  array of (3,n,m) shape, but need (n,m,3)
#     c = c.swapaxes(0,2)
#     return c


if __name__ == "__main__":

    h5file = "vx_id16a_A.h5"

    f = h5py.File(h5file,'r')

    arr1 = f["uptomode0/Wcomplex/image_data"].value
    x = f["uptomode0/Wcomplex/axis_x"].value
    y = f["uptomode0/Wcomplex/axis_y"].value
    f.close()
    #
    # print(arr1.shape,x.shape,y.shape)
    #
    #
    # plot_image(numpy.angle(arr1.T),x,y,cmap='hsv',
    #             xrange=[-0.05,0.05],xtitle="X [mm]",
    #             yrange=[-0.05,0.05],ytitle="Y [mm]",)


    # N = 1024
    # x, y = np.ogrid[-4:4:N*1j, -4:4:N*1j]
    # z = x + 1j*y

    # img = Complex2HSV(z, 0, 4)
    # img = Complex2HSV(arr1, 10,50)

    # Plot the array "A" using colorize

    # cmap = plt.cm.RdYlBu
    cmap = plt.cm.hsv
    phases = numpy.angle(arr1)


    colors = Normalize(phases.min(),phases.max(),clip=True)(phases)
    print(">>>> colors 0",colors.shape)
    colors = cmap(colors)
    print(">>>> colors 1",colors.shape)




    weights = numpy.abs(arr1)**2
    print("Extrema for weights: ",weights.min(),weights.max())


    rmin = 1e23
    rmax = weights.max()
    weights = np.where(weights < rmin, rmin, weights)
    weights = np.where(weights > rmax, rmax, weights)


    print("Extrema for weights: ",weights.min(),weights.max())
    weights = numpy.log10(weights)
    print("Extrema for weights: ",weights.min(),weights.max())

    weights -= weights.min()
    weights /= weights.max()

    print("Extrema for weights: ",weights.min(),weights.max())

    # rmin = 0.5
    # rmax = 0.6
    # weights = np.where(weights < rmin, rmin, weights)
    # weights = np.where(weights > rmax, rmax, weights)
    # weights -= weights.min()
    # weights /= weights.max()




    print(">>>",colors.shape,phases.shape,weights.shape)


    colors[..., -1] = weights
    print(">>>> colors 2",colors.shape)


    import pylab as plt
    plt.imshow(colors, interpolation='none',cmap=cmap)
    plt.colorbar()
    plt.show()

    # plot_image(weights.T)
