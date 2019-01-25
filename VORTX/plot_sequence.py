import numpy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize

import h5py


def plot_with_transparency_nine(arr1_list,DISTANCES,extent=(-75,75,-15,15),delta=6,
                                set_xlim=None,set_ylim=None,point_coordinates=None,
                                show=True,savefig=None):

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


    fig=plt.figure(figsize=(8, 8))
    fig.subplots_adjust(hspace=0.01, wspace=0.01)
    columns = int(numpy.sqrt(DISTANCES.size))
    rows = columns
    ax_list = []
    for i in range(1, columns*rows +1):
        img = colors_list[i-1]
        ax = fig.add_subplot(rows, columns, i)
        ax_list.append(ax)
        ax.imshow(img, interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')

        if i != 1 and i != 4 and i != 7:
            ax.yaxis.set_major_formatter(plt.NullFormatter())
        else:
            plt.ylabel("Y [$\mu$m]")

        if i  < 7:
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        else:
            plt.xlabel("X [$\mu$m]")

        ax.text(-70, -70, "$D$=%5.2f m"%DISTANCES[i-1])#, bbox={'facecolor': 'white', 'pad': 10})

    for i in range(len(ax_list)):
        ax_list[i].set_xlim( set_xlim )
        ax_list[i].set_ylim( set_ylim )

    if point_coordinates is not None:
        for i in range(len(ax_list)):
            ax_list[i].plot(point_coordinates[0],point_coordinates[1],marker='o',mfc='none')


    # plt.show()


    # fig, ax = plt.subplots(5, 5, figsize=(8, 8))
    # print(ax[0])
    # fig.subplots_adjust(hspace=0, wspace=0)
    #
    #
    # ax[0,0].imshow(colors_list[0], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    # ax[0,0].xaxis.set_major_formatter(plt.NullFormatter())
    #
    # ax[1].imshow(colors_list[1], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    # ax[1].xaxis.set_major_formatter(plt.NullFormatter())
    #
    # ax[2].imshow(colors_list[2], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')
    # ax[2].xaxis.set_major_formatter(plt.NullFormatter())
    #
    # ax[3].imshow(colors_list[3], interpolation='none',cmap=cmap,aspect='equal',extent=extent, origin='lower')



    if savefig is not None:
        plt.savefig(savefig)
        print("File written to disk ",savefig)

    if show: plt.show()


if __name__ == "__main__""":

    up_to_mode = 19
    h5file_root = "vx_id16a_C5_propagated_neighbour_mode%04d"%up_to_mode
    h5file = h5file_root+".h5"

    f = h5py.File(h5file,'r')

    DISTANCES = f["DISTANCES"].value
    list_with_images = []
    for i in range(DISTANCES.size):
        print("trying to read: ", "uptomode%04d_%04d/Wcomplex/image_data"%(up_to_mode,i))
        list_with_images.append(f["uptomode%04d_%04d/Wcomplex/image_data"%(up_to_mode,i)].value)

        print(list_with_images[i].shape)

    x = f["uptomode%04d_0000/Wcomplex/axis_x"%up_to_mode].value * 1e3 # in microns
    y = f["uptomode%04d_0000/Wcomplex/axis_y"%up_to_mode].value * 1e3 # in microns
    point_coordinates = f["r2"].value

    f.close()


    # from srxraylib.plot.gol import plot_image
    # plot_image(numpy.angle(list_with_images[1].T),x,y)

    # plot_with_transparency_four(list_with_images,DISTANCES,extent=(-25,25,-10,10),delta=8,savefig=h5file_root+".png")
    plot_with_transparency_nine(list_with_images,DISTANCES,extent=(x[0],x[-1],y[0],y[-1]),delta=8,
                point_coordinates=point_coordinates*1e6,set_xlim=[-75,75],set_ylim=None,
                savefig="/tmp/"+h5file_root+".png")

    #
    #
    # w=10
    # h=10
    # fig=plt.figure(figsize=(8, 8))
    # columns = 4
    # rows = 5
    # for i in range(1, columns*rows +1):
    #     img = np.random.randint(10, size=(h,w))
    #     fig.add_subplot(rows, columns, i)
    #     plt.imshow(img)
    # plt.show()