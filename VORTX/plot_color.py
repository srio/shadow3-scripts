import numpy
import numpy as np
import h5py
from srxraylib.plot.gol import plot_image,plot_surface

# see https://stackoverflow.com/questions/17044052/mathplotlib-imshow-complex-2d-array

from colorsys import hls_to_rgb
from matplotlib.colors import hsv_to_rgb

import pylab as plt
from matplotlib.colors import Normalize


        #arr1,extent=(-75,75,-15,15),delta=6,cmap=None,show=True):
def plot_with_transparency_one(arr0,extent=None,delta=6,cmap=None,show=True):

    from colorsys import hls_to_rgb
    from matplotlib.colors import hsv_to_rgb



    arr1 = arr0.T

    cmap = plt.cm.hsv

    phases = numpy.angle(arr1)


    colors = Normalize(phases.min(),phases.max(),clip=True)(phases)
    colors = cmap(colors)


    weights = numpy.abs(arr1)**2


    rmax = weights.max()
    rmin = rmax/(10**delta) # 1e23
    weights = numpy.where(weights < rmin, rmin, weights)
    weights = numpy.where(weights > rmax, rmax, weights)

    weights = numpy.log10(weights)

    weights -= weights.min()
    weights /= weights.max()

    colors[..., -1] = weights

    fig = plt.figure()

    plt.imshow(colors, interpolation='none',cmap=cmap,aspect='equal',origin='lower')


    plt.ylim( (x.min(),x.max()) )
    plt.ylim( (y.min(),y.max()) )

    # if cmap is not None:
    #     plt.colorbar()

    if show: plt.show()



def plot_with_transparency_four(arr1_list,x,y,extent=(-75,75,-15,15),delta=6,show=True,savefig="/tmp/1.png"):

    cmap = plt.cm.hsv


    colors_list = []

    for i,arr1 in enumerate(arr1_list):
        phases = numpy.angle(arr1)


        colors = Normalize(phases.min(),phases.max(),clip=True)(phases)
        colors = cmap(colors)


        weights = numpy.abs(arr1)**2
        print("Extrema for weights: ",weights.min(),weights.max())


        rmax = weights.max()
        rmin = rmax/(10**delta) # 1e23
        weights = np.where(weights < rmin, rmin, weights)
        weights = np.where(weights > rmax, rmax, weights)

        weights = numpy.log10(weights)

        weights -= weights.min()
        weights /= weights.max()


        colors[..., -1] = weights

        colors_list.append(colors)


    fig, ax = plt.subplots(4, 1, figsize=(8, 8))
    print(ax[0])
    fig.subplots_adjust(hspace=0, wspace=0)


    ax[0].imshow(colors_list[0], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[0].xaxis.set_major_formatter(plt.NullFormatter())

    ax[1].imshow(colors_list[1], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[1].xaxis.set_major_formatter(plt.NullFormatter())

    ax[2].imshow(colors_list[2], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[2].xaxis.set_major_formatter(plt.NullFormatter())

    ax[3].imshow(colors_list[3], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')


    plt.xlim( (x.min(),x.max()) )
    plt.ylim( (y.min(),y.max()) )

    plt.xlabel("X [$\mu$m]")
    plt.ylabel("Y [$\mu$m]")


    if savefig is not None:
        plt.savefig(savefig)
        print("File written to disk ",savefig)

    if show: plt.show()


if __name__ == "__main__":

    h5file_root = "vx_id16a_C5_propagated"
    h5file = h5file_root+".h5"

    f = h5py.File(h5file,'r')

    arr1 = f["uptomode0000/Wcomplex/image_data"].value
    arr2 = f["uptomode0009/Wcomplex/image_data"].value
    arr3 = f["uptomode0099/Wcomplex/image_data"].value
    arr4 = f["uptomode0999/Wcomplex/image_data"].value

    x    = f["uptomode0000/Wcomplex/axis_x"].value
    y    = f["uptomode0000/Wcomplex/axis_y"].value
    f.close()

    print("Data limits X: ",x.min(),x.max()," Y: ",y.min(),y.max())
    plot_with_transparency_four([arr1,arr2,arr3,arr4],x,y,delta=8,savefig=h5file_root+".png")#,extent=(-25,25,-10,10),

    # plot_with_transparency_one(arr2,cmap='hsv',extent=(-25,25,-10,10))

    # tmp = numpy.angle(arr3)
    # print(tmp.shape)
    # plot_image(tmp.T,cmap='hsv')

    #
    # cmap = plt.cm.hsv
    # phases = numpy.angle(arr1)
    #
    #
    # colors = Normalize(phases.min(),phases.max(),clip=True)(phases)
    # # colors = (phases)
    # print(">>>> colors 0",colors.shape)
    # colors = cmap(colors)
    # print(">>>> colors 1",colors.shape)
    #
    #
    #
    # delta = 6
    # weights = numpy.abs(arr1)**2
    # print("Extrema for weights: ",weights.min(),weights.max())
    #
    #
    # rmax = weights.max()
    # rmin = rmax/(10**delta) # 1e23
    # weights = np.where(weights < rmin, rmin, weights)
    # weights = np.where(weights > rmax, rmax, weights)
    #
    #
    # print("Extrema for weights: ",weights.min(),weights.max())
    # weights = numpy.log10(weights)
    # print("Extrema for weights: ",weights.min(),weights.max())
    #
    # weights -= weights.min()
    # weights /= weights.max()
    #
    # # weights = weights * 2 * np.pi + np.pi
    #
    # print("Extrema for weights: ",weights.min(),weights.max())
    #
    # # rmin = 0.5
    # # rmax = 0.6
    # # weights = np.where(weights < rmin, rmin, weights)
    # # weights = np.where(weights > rmax, rmax, weights)
    # # weights -= weights.min()
    # # weights /= weights.max()
    #
    #
    #
    #
    # print(">>>",colors.shape,phases.shape,weights.shape)
    #
    #
    # colors[..., -1] = weights
    # print(">>>> colors 2",colors.shape)
    #
    # fig = plt.figure()
    # plt.subplot(411)
    #
    # plt.imshow(colors, interpolation='none',cmap=cmap,aspect='equal',
    #            extent=(-75,75,-15,15))
    #
    #
    # # plt.colorbar()
    # # plt.imshow(weights, interpolation='none',cmap=cmap)
    # # plt.ylim( (x.min(),x.max()) )
    # # plt.ylim( (y.min(),y.max()) )
    #
    # plt.show()
    #
    # # plot_image(colors.T,aspect='auto')
