import numpy
import numpy as np
import h5py
from srxraylib.plot.gol import plot_image,plot_surface,plot

# see https://stackoverflow.com/questions/17044052/mathplotlib-imshow-complex-2d-array

from colorsys import hls_to_rgb
from matplotlib.colors import hsv_to_rgb

import pylab as plt
from matplotlib.colors import Normalize
import scipy.constants as codata

        #arr1,extent=(-75,75,-15,15),delta=6,cmap=None,show=True):
def plot_with_transparency_one(arr0,x,y,extent=None,delta=6,cmap=plt.cm.hsv,show=True):

    from colorsys import hls_to_rgb
    from matplotlib.colors import hsv_to_rgb

    print(x)
    extent=(x[0],x[-1],y[0],y[-1])

    arr1 = arr0 #.T

    # cmap = plt.cm.grey
    # cmap = plt.cm.hsv

    # phases = numpy.angle(arr1)

    X = numpy.outer(x,numpy.ones_like(y)).T * 1e-6 # from um to m
    Y = numpy.outer(numpy.ones_like(x),y).T * 1e-6 # from um to m
    energy = 17226.0
    wavelength = codata.h * codata.c / codata.e / energy
    background = numpy.exp(2j * numpy.pi / wavelength * (X * X + Y * Y)  / 2 / 5.0)

    phases = numpy.angle(arr1 * background)

    plot(x, numpy.angle(arr1)[y.size//2,:],
         x, numpy.angle(background)[y.size//2,:],
         x, phases[y.size//2,:],legend=["comsyl","background","diff"],xtitle="x",show=False)

    plot(y, numpy.angle(arr1)[:,x.size//2],
         y, numpy.angle(background)[:,x.size//2],
         y, phases[:,x.size//2],legend=["comsyl","background","diff"],xtitle="y")

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

    plt.imshow(colors, interpolation='none',cmap=cmap,aspect='equal',origin='lower',extent=extent)


    plt.ylim( (x.min(),x.max()) )
    plt.ylim( (y.min(),y.max()) )

    # if cmap is not None:
    #     plt.colorbar()

    if show: plt.show()



def plot_with_transparency_four(arr1_list,x,y,delta=6,show=True,savefig="/tmp/1.png",point_coordinates=None,
                                propagation_distance=0,set_xlim=None,set_ylim=None,
                                correct_background=False,cmap=plt.cm.hsv,
                                plot_intensity=False):

    extent=(x[0],x[-1],y[0],y[-1])

    # cmap = plt.cm.hsv


    colors_list = []

    if correct_background:
        if propagation_distance == 0:
            background = 1.0
        else:
            X = numpy.outer(x, numpy.ones_like(y)).T * 1e-6  # from um to m
            Y = numpy.outer(numpy.ones_like(x), y).T * 1e-6  # from um to m
            energy = 17226.0
            wavelength = codata.h * codata.c / codata.e / energy
            background = numpy.exp(2j * numpy.pi / wavelength * (X * X + Y * Y) / 2 / propagation_distance)
    else:
        background = 1.0


    for i,arr1 in enumerate(arr1_list):


        phases = numpy.angle(arr1*background)

        weights = numpy.abs(arr1)**2
        print("Extrema for weights: ",weights.min(),weights.max())

        if plot_intensity:
            w1 = weights.copy() * 0 + 1
            colors = Normalize(0, 1, clip=True)(w1)
        else: # phase, default
            colors = Normalize(phases.min(),phases.max(),clip=True)(phases)



        colors = cmap(colors)


        rmax = weights.max()
        rmin = rmax/(10**delta) # 1e23
        weights = np.where(weights < rmin, rmin, weights)
        weights = np.where(weights > rmax, rmax, weights)

        weights = numpy.log10(weights)

        weights -= weights.min()
        weights /= weights.max()


        colors[..., -1] = weights #* 0 + 1

        colors_list.append(colors)


    # fig, ax = plt.subplots(4, 1, figsize=(4, 8))
    if propagation_distance == 0:
        fig, ax = plt.subplots(4, 1, figsize=(2.5, 9))
    else: # if propagation_distance == 30:
        fig, ax = plt.subplots(4, 1, figsize=(2.5, 8.9))

    print(ax[0])
    fig.subplots_adjust(hspace=0, wspace=0)
    # fig.tight_layout()


    ax[0].imshow(colors_list[0], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[0].xaxis.set_major_formatter(plt.NullFormatter())
    ax[0].yaxis.set_major_formatter(plt.NullFormatter())


    ax[1].imshow(colors_list[1], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[1].xaxis.set_major_formatter(plt.NullFormatter())
    ax[1].yaxis.set_major_formatter(plt.NullFormatter())

    ax[2].imshow(colors_list[2], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[2].xaxis.set_major_formatter(plt.NullFormatter())
    ax[2].yaxis.set_major_formatter(plt.NullFormatter())

    ax[3].imshow(colors_list[3], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    ax[3].yaxis.set_major_formatter(plt.NullFormatter())


    for i in range(4):
        ax[i].set_xlim( set_xlim )
        ax[i].set_ylim( set_ylim )


    if propagation_distance == 0:
        if point_coordinates is not None:
            for i in range(4):
                ax[i].plot(point_coordinates[0],point_coordinates[1],marker='o',mfc='none')
                ax[i].tick_params(labelsize=12)
                ax[i].yaxis.set_ticks([-30,0,30])
                ax[i].xaxis.set_ticks([-30, 0, 30])
    else: # if propagation_distance == 30:
        if point_coordinates is not None:
            for i in range(4):
                ax[i].plot(point_coordinates[0],point_coordinates[1],marker='o',mfc='none')
                ax[i].tick_params(labelsize=12)
                # ax[i].yaxis.set_ticks([-30,0,30])
                # ax[i].xaxis.set_ticks([-30, 0, 30])

    plt.xlabel("X [$\mu$m]")
    # plt.ylabel("Y [$\mu$m]")

    # plt.xticks([-25,0,25])
    # plt.yticks([-25,0,25])


    if savefig is not None:
        plt.savefig(savefig,dpi=600)
        print("File written to disk ",savefig)

    if show: plt.show()


if __name__ == "__main__":

    remove_spherical_phase = True
    cmap = plt.cm.hsv # YlGnBu #gray #hsv
    # cmap = plt.cm.binary

    # h5file_root = "vx_id16a_A"
    # propagation_distance = 0

    # h5file_root = "vx_id16a_B"
    # propagation_distance = 0
    #
    h5file_root = "vx_id16a_C"
    propagation_distance = 0

    # h5file_root = "vx_id16a_C1_propagated"
    # propagation_distance = 1

    # h5file_root = "vx_id16a_C5_propagated"
    # propagation_distance = 5

    # h5file_root = "vx_id16a_C_propagated"  # note missing "30"
    # propagation_distance = 30


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


    if propagation_distance == 30:
        plot_with_transparency_four([arr1,arr2,arr3,arr4], x, y, delta=8, savefig="/tmp/"+h5file_root+".png",
                                    propagation_distance=propagation_distance,
                                    point_coordinates=point_coordinates, set_xlim=[-400,400], set_ylim=[-400,400],
                                    correct_background=remove_spherical_phase,cmap=cmap,show=True)
        # plot_with_transparency_one(arr1,x,y,delta=8)
    elif propagation_distance == 5:
        plot_with_transparency_four([arr1,arr2,arr3,arr4], x, y, delta=8, savefig="/tmp/"+h5file_root+".png",
                                    propagation_distance=propagation_distance,
                                    point_coordinates=point_coordinates, set_xlim=[-70,70], set_ylim=[-70,70],
                                    correct_background=remove_spherical_phase,cmap=cmap,show=True,)

        # plot_with_transparency_one(arr1, x, y, delta=8)
    elif propagation_distance == 1:
        plot_with_transparency_four([arr1,arr2,arr3,arr4], x, y, delta=8, savefig="/tmp/"+h5file_root+".png",
                                    propagation_distance=propagation_distance,
                                    point_coordinates=point_coordinates, set_xlim=[-35,35], set_ylim=[-35,35],
                                    correct_background=remove_spherical_phase,cmap=cmap,show=True)
    elif propagation_distance == 0:
        plot_with_transparency_four([arr1, arr2, arr3, arr4], x, y, delta=8, savefig="/tmp/" + h5file_root + ".png",
                                    propagation_distance=propagation_distance,
                                    point_coordinates=point_coordinates, set_xlim=[-50, 50],
                                    correct_background=remove_spherical_phase,cmap=cmap,
                                    plot_intensity=False,
                                    # set_xlim=[-575,575],set_ylim=[-575,575]
                                    )


    else:
        raise("Not implemented")


