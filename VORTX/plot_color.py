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
# def plot_with_transparency_one(arr0,extent=None,delta=6,cmap=None,show=True):
#
#     from colorsys import hls_to_rgb
#     from matplotlib.colors import hsv_to_rgb
#
#
#
#     arr1 = arr0.T
#
#     cmap = plt.cm.hsv
#
#     phases = numpy.angle(arr1)
#
#
#     colors = Normalize(phases.min(),phases.max(),clip=True)(phases)
#     colors = cmap(colors)
#
#
#     weights = numpy.abs(arr1)**2
#
#
#     rmax = weights.max()
#     rmin = rmax/(10**delta) # 1e23
#     weights = numpy.where(weights < rmin, rmin, weights)
#     weights = numpy.where(weights > rmax, rmax, weights)
#
#     weights = numpy.log10(weights)
#
#     weights -= weights.min()
#     weights /= weights.max()
#
#     colors[..., -1] = weights
#
#     fig = plt.figure()
#
#     plt.imshow(colors, interpolation='none',cmap=cmap,aspect='equal',origin='lower')
#
#
#     plt.ylim( (x.min(),x.max()) )
#     plt.ylim( (y.min(),y.max()) )
#
#     # if cmap is not None:
#     #     plt.colorbar()
#
#     if show: plt.show()



def plot_with_transparency_four(arr1_list,x,y,delta=6,show=True,savefig="/tmp/1.png",point_coordinates=None,
                                set_xlim=None,set_ylim=None):

    extent=(x[0],x[-1],y[0],y[-1])

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


    fig, ax = plt.subplots(4, 1, figsize=(4, 8))
    print(ax[0])
    fig.subplots_adjust(hspace=0, wspace=0)


    ax[0].imshow(colors_list[0], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[0].xaxis.set_major_formatter(plt.NullFormatter())


    ax[1].imshow(colors_list[1], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[1].xaxis.set_major_formatter(plt.NullFormatter())

    ax[2].imshow(colors_list[2], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[2].xaxis.set_major_formatter(plt.NullFormatter())

    ax[3].imshow(colors_list[3], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')

    for i in range(4):
        ax[i].set_xlim( set_xlim )
        ax[i].set_ylim( set_ylim )

    if point_coordinates is not None:
        for i in range(4):
            ax[i].plot(point_coordinates[0],point_coordinates[1],marker='o',mfc='none')

    plt.xlabel("X [$\mu$m]")
    plt.ylabel("Y [$\mu$m]")


    if savefig is not None:
        plt.savefig(savefig)
        print("File written to disk ",savefig)

    if show: plt.show()


if __name__ == "__main__":

    h5file_root = "vx_id16a_C_propagated"



    #
    #
    #

    h5file = h5file_root+".h5"

    f = h5py.File(h5file,'r')

    point_coordinates = f["r2"].value
    arr1 = f["uptomode0000/Wcomplex/image_data"].value
    arr2 = f["uptomode0009/Wcomplex/image_data"].value
    arr3 = f["uptomode0099/Wcomplex/image_data"].value
    arr4 = f["uptomode0999/Wcomplex/image_data"].value

    x    = f["uptomode0000/Wcomplex/axis_x"].value * 1e3 # in microns
    y    = f["uptomode0000/Wcomplex/axis_y"].value * 1e3 # in microns
    f.close()

    print("Data limits [um] X: ",x.min(),x.max()," Y: ",y.min(),y.max())
    print("Point coordinates [um] X: ",point_coordinates[0]*1e6," Y: ",point_coordinates[1]*1e6)
    point_coordinates *= 1e6
    # if ( (numpy.abs(point_coordinates[0]) > extent[1]) or (numpy.abs(point_coordinates[1]) > extent[3])):
    #     point_coordinates_plot = None
    # else:
    # #     point_coordinates_plot = point_coordinates
    # extent=(-25,25,-15,15)
    # extent=(x[0],x[-1],y[0],y[-1])
    # point_coordinates = None
    plot_with_transparency_four([arr1,arr2,arr3,arr4],x,y,delta=8,savefig=h5file_root+".png",
                                point_coordinates=point_coordinates) #,set_xlim=[-575,575],set_ylim=[-575,575])


