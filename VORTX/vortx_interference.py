from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter

from srxraylib.plot.gol import plot_image, plot

import h5py

from vortx_propagate import AFpropagated #, W_at_x2x2, propagate, apply_two_apertures
#
# from plot_color import plot_with_transparency_one


import pylab as plt
from matplotlib.colors import Normalize, ListedColormap
import matplotlib.patches as patches

########################################################################################################################



def plot_image_with_transparency(
            *positional_parameters,title="TITLE",xtitle=r"X",ytitle=r"Y",
            transparency_log=True,delta=6,
            xrange=None, yrange=None,cmap=None,aspect=None,add_colorbar=True,figsize=None,show=True,
            patch_shape=None,
            patch1_center=None,patch1_width=None,
            patch2_center=None,patch2_width=None ):

    n_arguments = len(positional_parameters)
    if n_arguments == 1:
        raise Exception("At least two arguments required!")
    elif n_arguments == 2:
        z = positional_parameters[0].T
        weights = positional_parameters[1].T
        x = numpy.arange(0,z.shape[1])
        y = numpy.arange(0,z.shape[0])
    elif n_arguments == 3:
        z = positional_parameters[0].T
        weights = positional_parameters[1].T
        x = positional_parameters[2]
        y = numpy.arange(0,z.shape[0])
    elif n_arguments == 4:
        z = positional_parameters[0].T
        weights = positional_parameters[1].T
        x = positional_parameters[2]
        y = positional_parameters[3]
    else:
        raise Exception("Bad number of inputs")

    extent = [x.min(),x.max(),y.min(),y.max()]
    if xrange is not None:
        extent[0] = xrange[0]
        extent[1] = xrange[1]
    if yrange is not None:
        extent[2] = yrange[0]
        extent[3] = yrange[1]



    fig = plt.figure(figsize=figsize)


    if isinstance(cmap,ListedColormap):
        cmap1 = cmap
    elif isinstance(cmap,str):
        cmap1 = plt.cm.get_cmap(name=cmap, lut=None) #plt.cm.hsv
    else:
        cmap1 = plt.cm.get_cmap(name=None, lut=None)


    # colors = Normalize(z.min(),z.max(),clip=True)(z)
    colors = Normalize(vmin=-numpy.pi, vmax=numpy.pi,clip=True)(z)


    print(">>>>1",colors.shape)
    colors = cmap1(colors)
    print(">>>>2",colors.shape)

    # weights = zt # numpy.abs(arr1)**2

    rmax = weights.max()
    rmin = rmax/(10**delta) # 1e23
    weights = numpy.where(weights < rmin, rmin, weights)
    weights = numpy.where(weights > rmax, rmax, weights)

    if transparency_log:
        weights = numpy.log10(weights)

    weights -= weights.min()
    weights /= weights.max()

    colors[..., -1] = weights

    # fig = plt.figure()

    plt.imshow(colors, interpolation='none',origin='lower',extent=[x[0],x[-1],y[0],y[-1]],cmap=cmap1,aspect=aspect)

    ax = fig.gca()
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)



    if patch_shape is not None:

        if patch_shape == "Rectangle":
            rect1 = patches.Rectangle(
                    (patch1_center[0]-0.5*patch1_width[0],patch1_center[1]-0.5*patch1_width[1]),
                        patch1_width[0],patch1_width[1],
                                     linewidth=1,edgecolor='r',facecolor='none')
            ax.add_patch(rect1)


            if patch2_center is not None:
                rect2 = patches.Rectangle(
                        (patch2_center[0]-0.5*patch2_width[0],patch2_center[1]-0.5*patch2_width[1]),
                            patch2_width[0],patch2_width[1],
                                         linewidth=1,edgecolor='k',facecolor='none')

                ax.add_patch(rect2)

        elif patch_shape == "Ellipse":

            ell1 = patches.Ellipse(
                    (patch1_center[0],patch1_center[1]),
                        patch1_width[0],patch1_width[1],
                                     linewidth=1,edgecolor='r',facecolor='none')
            ax.add_patch(ell1)

            if patch2_center is not None:
                ell2 = patches.Ellipse(
                        (patch2_center[0],patch2_center[1]),
                            patch2_width[0],patch2_width[1],
                                         linewidth=1,edgecolor='k',facecolor='none')
                ax.add_patch(ell2)

    plt.title(title)

    plt.xlim( xrange )
    plt.ylim( yrange )

    #
    # if add_colorbar:
    #     plt.colorbar()

    if show:
        plt.show()

def plot_color_table(orientation='horizontal'):

    import matplotlib as mpl
    if orientation == 'horizontal':
        fig = plt.figure(figsize=(10,2))
        ax1 = fig.add_axes([0.05, 0.40, 0.9, 0.25])
    else:
        raise NotImplementedError


    # Make a figure and axes with dimensions as desired.
    # fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_axes([0.05, 0.40, 0.9, 0.25])
    # ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])
    # ax3 = fig.add_axes([0.05, 0.15, 0.9, 0.15])

    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    cmap = mpl.cm.hsv
    norm = mpl.colors.Normalize(vmin=-numpy.pi, vmax=numpy.pi)

    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                    norm=norm,
                                    orientation=orientation)
    cb1.set_label('Phase [rad]')
    #
    #
    # ################################

    plt.show()

if __name__ == "__main__":

    write_h5 = True
    filename_ebs="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"

    point = "J" #"C5"

    distance = 30.0
    index_max = 19
    zoom = (6.0,16.0)

    if point == "A":
        coordinate_x = 0.0
        coordinate_y = 0.0
    elif point == "B":
        coordinate_x = 57.0e-6
        coordinate_y = 148.2e-6
    elif point == "C":
        coordinate_x = 125.47e-6
        coordinate_y = 312.78e-6
    elif point == "D":
        coordinate_x = -7e-6
        coordinate_y = 31e-6
    elif point == "C1":
        distance = 1.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        zoom = (1.0,1.0)
    elif point == "C5":
        distance = 5.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        zoom = (2.0,2.0)
    elif point == "C100":
        distance = 100.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        zoom = (10.0,40.0)
    elif point == "J":
        distance = 5.0
        coordinate_x = -10.0e-6 # 0.33 * 125.47e-6 / 30.0 * distance
        coordinate_y = -25.0e-6 # 0.33 * 312.78e-6 / 30.0 * distance
        zoom = (2.0,2.0)
    else:
        raise Exception("Point not found!")


    #
    # slits
    #
    patch_shape = "Ellipse"
    center1 = [coordinate_x,coordinate_y]
    width1 = [4e-6,4e-6]
    center2 = [10e-6,25e-6] # [5e-6,25e-6]
    width2 = width1

    patch1_center = [1e6*center1[0],1e6*center1[1]]
    patch1_width  = [ 1e6*width1[0], 1e6*width1[1]]

    patch2_center = [1e6*center2[0],1e6*center2[1]]
    patch2_width  = [ 1e6*width2[0], 1e6*width2[1]]




    #
    # load CSD
    #


    af  = CompactAFReader.initialize_from_file(filename_ebs)


    #
    # get indices
    #

    # first propagate a few modes only to check there are no errors
    afp = AFpropagated.propagate(af,distance=distance,index_max=1,zoom=zoom)

    if write_h5:
        h5file = "vx_id16a_%s_propagated.h5"%point

    print("X: start, step, points",afp.x_coordinates()[0],afp.x_coordinates()[1] - afp.x_coordinates()[0],afp.x_coordinates().size)
    print("Y: start, step, points",afp.y_coordinates()[0],afp.y_coordinates()[1] - afp.y_coordinates()[0],afp.y_coordinates().size)
    step_x = afp.x_coordinates()[1] - afp.x_coordinates()[0]
    step_y = afp.y_coordinates()[1] - afp.y_coordinates()[0]
    origin_x = afp.x_coordinates()[0]
    origin_y = afp.y_coordinates()[0]

    index_x2 = int((coordinate_x - origin_x) / step_x)
    index_y2 = int((coordinate_y - origin_y) / step_y)


    print("Using indices: ",index_x2,index_y2," out of ",afp.x_coordinates().size,afp.y_coordinates().size,
          "ratio: ",index_x2/afp.x_coordinates().size,index_y2/afp.y_coordinates().size)


    if write_h5:
        h5w = AFpropagated.h5_initialize("tmp.h5")
    #
    # propagate
    #

    # now propagate all modes
    afp = AFpropagated.propagate(af,distance=distance,index_max=index_max,zoom=zoom)


    # plot CSD with slits
    tmp = afp.W_at_x2x2(index_x2=index_x2,index_y2=index_y2,index_max=index_max)
    x = afp.x_coordinates()
    y = afp.y_coordinates()

    plot_image_with_transparency(numpy.angle(tmp),numpy.abs(tmp)**2,1e6*x,1e6*y,
                title="uptomode%04d"%index_max,
                xtitle="X [um, %d pixels]"%x.size,ytitle="Y [um, %d pixels]"%y.size,cmap='hsv',show=False,
                xrange=[-150,150],yrange=[-100,100],aspect='equal',add_colorbar=False,patch_shape=patch_shape,
                patch1_center=patch1_center,patch1_width=patch1_width,
                patch2_center=patch2_center,patch2_width=patch1_width)

    if write_h5:
        afp.h5w = h5w
        afp.h5_W_at_x2x2(index_x2=index_x2,index_y2=index_y2,index_max=index_max)

    #
    # slits
    #

    afp_cut = afp.apply_two_apertures(index_max=index_max,patch_shape=patch_shape,
                                         center1=center1,width1=width1,center2=center2,width2=width2)

    tmp = afp_cut.W_at_x2x2(index_x2=index_x2,index_y2=index_y2,index_max=index_max)
    x = afp.x_coordinates()
    y = afp.y_coordinates()

    print("Total intensiy before slits: ",afp.get_intensity().sum())
    print("Total intensiy after slits: ",afp_cut.get_intensity().sum())



    # plot CSD with slits
    tmp = afp_cut.W_at_x2x2(index_x2=index_x2,index_y2=index_y2,index_max=index_max)
    x = afp.x_coordinates()
    y = afp.y_coordinates()

    plot_image_with_transparency(numpy.angle(tmp),numpy.abs(tmp)**2,1e6*x,1e6*y,
                title="uptomode%04d"%index_max,
                xtitle="X [um, %d pixels]"%x.size,ytitle="Y [um, %d pixels]"%y.size,cmap='hsv',show=False,
                xrange=[-150,150],yrange=[-100,100],aspect='equal',add_colorbar=False,
                patch1_center=patch1_center,patch1_width=patch1_width,
                patch2_center=patch2_center,patch2_width=patch1_width)

    #
    # propagate again
    #

    # afpp = propagate(afp_cut,distance=15,index_max=index_max,zoom=(2.0,2.0))
    afpp = AFpropagated.propagate(afp_cut,distance=35,index_max=index_max,zoom=(2.0,5.0))


    tmp = afpp.get_intensity(index_max) #W_at_x2x2(afpp,index_x2=index_x2,index_y2=index_y2,index_max=index_max)
    x = afpp.x_coordinates()
    y = afpp.y_coordinates()


    plot_image(numpy.array(tmp),1e6*x,1e6*y,title="Two slits interference",#xrange=[-150,150],yrange=[-100,100],
                xtitle="X [um, %d pixels]"%x.size,ytitle="Y [um, %d pixels]"%y.size,cmap='hsv',show=False,
                aspect='equal')

    plot_image(numpy.log10(tmp),1e6*x,1e6*y,title="Two slits interference - LOG!!",#xrange=[-150,150],yrange=[-100,100],
                xtitle="X [um, %d pixels]"%x.size,ytitle="Y [um, %d pixels]"%y.size,cmap='hsv',show=False,
                aspect='equal')
    #
    #




    if write_h5:

        afpp.h5w = h5w

        afpp.h5_get_intensity(index_max,"intensity_at_image")

        for i in range(afpp.number_modes()):
            print("adding to h5 file mode : ",i)
            afpp.h5_add_mode(i)


    plt.show()


    print("coordinates [um]",1e6*coordinate_x,1e6*coordinate_y)