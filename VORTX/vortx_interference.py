from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter

from srxraylib.plot.gol import plot_image, plot

import h5py

from vortx_propagate import AFpropagated, W_at_x2x2, propagate
#
# from plot_color import plot_with_transparency_one


import pylab as plt
from matplotlib.colors import Normalize, ListedColormap
########################################################################################################################



def plot_image_with_transparency(
            *positional_parameters,title="TITLE",xtitle=r"X",ytitle=r"Y",
            transparency_log=True,delta=6,
            xrange=None, yrange=None,cmap=None,aspect=None,add_colorbar=True,figsize=None,show=True):

    n_arguments = len(positional_parameters)
    if n_arguments == 1:
        raise Exception("At least two arguments required!")
    elif n_arguments == 2:
        z = positional_parameters[0].T
        zt = positional_parameters[1].T
        x = numpy.arange(0,z.shape[0])
        y = numpy.arange(0,z.shape[1])
    elif n_arguments == 3:
        z = positional_parameters[0].T
        zt = positional_parameters[1].T
        x = positional_parameters[2]
        y = numpy.arange(0,z.shape[1])
    elif n_arguments == 4:
        z = positional_parameters[0].T
        zt = positional_parameters[1].T
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
    #
    # # cmap = plt.cm.Greys
    # plt.imshow(z.T,origin='lower',extent=[x[0],x[-1],y[0],y[-1]],cmap=cmap,aspect=aspect)
    # if add_colorbar:
    #     plt.colorbar()
    # ax = fig.gca()
    # ax.set_xlabel(xtitle)
    # ax.set_ylabel(ytitle)
    #
    # plt.title(title)
    #
    # plt.xlim( xrange )
    # plt.ylim( yrange )
    #
    # if show:
    #     plt.show()
    #
    # return fig



    if isinstance(cmap,ListedColormap):
        cmap1 = cmap
    elif isinstance(cmap,str):
        cmap1 = plt.cm.get_cmap(name=cmap, lut=None) #plt.cm.hsv
    else:
        cmap1 = plt.cm.get_cmap(name=None, lut=None)

    phases = z # numpy.angle(arr1)

    colors = Normalize(phases.min(),phases.max(),clip=True)(phases)
    colors = cmap1(colors)

    weights = zt # numpy.abs(arr1)**2

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

    plt.title(title)

    plt.xlim( xrange )
    plt.ylim( yrange )


    if add_colorbar:
        plt.colorbar()

    if show:
        plt.show()

if __name__ == "__main__":


    filename_ebs="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"

    point = "C5"

    distance = 30.0
    index_max = 9
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
    else:
        raise Exception("Point not found!")



    af  = CompactAFReader.initialize_from_file(filename_ebs)


    #
    # get indices
    #

    # first propagate a few modes only to check there are no errors
    afp = propagate(af,distance=distance,index_max=1,zoom=zoom)

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



    #
    # propagate
    #

    # now propagate all modes
    afp = propagate(af,distance=distance,index_max=index_max,zoom=zoom)


    tmp = W_at_x2x2(afp,index_x2=index_x2,index_y2=index_y2,index_max=index_max)


    #extent=(-25,25,-10,10),delta=8,
    x = afp.x_coordinates()
    y = afp.y_coordinates()

    print(">>>>>>>>>>",tmp.shape,x.shape,afp.y.shape)

    plot_image(numpy.angle(tmp),x,y,
                title="uptomode%04d"%index_max,
                xtitle="X [um]",ytitle="Y [um]",cmap='hsv',show=False)

    plot_image_with_transparency(numpy.angle(tmp),numpy.abs(tmp)**2,1e6*x,1e6*y,
                title="uptomode%04d"%index_max,
                xtitle="X [um]",ytitle="Y [um]",cmap=None,show=True,
                xrange=[-100,100],yrange=[-30,30],aspect='equal')

    plt.show()
    # plot_with_transparency_one(tmp,extent=(-25,25,-10,10),delta=8,show=True)


    #

    # h5w = H5SimpleWriter.initialize_file(h5file,creator="vortx_propagate.py")
    # # h5_initialize(h5file,creator="vortx.py")
    # h5w.add_key("r2_indices",[index_x2,index_y2])
    # h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]])
    #
    #
    # x_ebs = numpy.arange(afp.number_modes())
    # y_ebs = numpy.abs(afp.occupation_array())
    #
    # # plot(x_ebs,y_ebs)
    #
    # if index_max > 100:
    #     t = numpy.array( (0,1,2,3,4,5,6,7,8,9,
    #                       19,29,39,49,59,69,79,89,99,
    #                       199,299,399,499,599,699,799,899,999,
    #                       1099) )
    # else:
    #     t = numpy.array( (0,1,2,3,4,5,6,7,8,9,
    #                       19,29,39,49,59,69,79,89,99, ))
    #
    #
    # # t = numpy.array( (0,1,2,3,4,5,6,7,8,9) )
    #
    # for index_max in t:
    #     print("Calculating mode index %d"%index_max)
    #     tmp = W_at_x2x2(afp,index_x2=index_x2,index_y2=index_y2,index_max=index_max)
    #
    #     h5w.create_entry("uptomode%04d"%index_max,nx_default="SpectralDensity")
    #     h5w.add_key("r2_indices",[index_x2,index_y2], entry_name="uptomode%04d"%index_max)
    #     h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]], entry_name="uptomode%04d"%index_max)
    #
    #     h5w.add_image(tmp,1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
    #                  entry_name="uptomode%04d"%index_max,
    #                  image_name="Wcomplex",title_x="X [mm]",title_y="Y [mm]")
    #
    #     h5w.add_image(numpy.absolute(tmp),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
    #                  entry_name="uptomode%04d"%index_max,
    #                  image_name="Wamplitude",title_x="X [mm]",title_y="Y [mm]")
    #
    #     h5w.add_image(numpy.angle(tmp),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
    #                  entry_name="uptomode%04d"%index_max,
    #                  image_name="Wphase",title_x="X [mm]",title_y="Y [mm]")
    #
    #     h5w.add_image(afp.get_intensity(index_max),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
    #                  entry_name="uptomode%04d"%index_max,
    #                  image_name="SpectralDensity",title_x="X [mm]",title_y="Y [mm]")
    #
    # print("File written to disk: ",h5file)