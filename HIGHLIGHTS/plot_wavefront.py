

import numpy

from srwlib import *
import sys
from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction
from comsyl.autocorrelation.AutocorrelationFunctionPropagator import AutocorrelationFunctionPropagator
from comsyl.parallel.utils import isMaster, barrier
from comsyl.utils.Logger import log


from comsyl.waveoptics.Wavefront import NumpyWavefront, SRWWavefront
from srxraylib.plot.gol import plot_image


def plot_filename(filename="",title="",show=True):

    #
    # wavefronts
    #


    # filename = "/users/srio/Working/paper-hierarchical/CODE-COMSYL/tmp/tmp0_chac_in_mode0.npz"
    # filename = "tmp/tmp0_hib3-3302_out.npz"
    a = NumpyWavefront.load(filename)

    print(a)

    # b = SRWWavefront(a)
    #
    # print(b)

    # print(a._e_field.shape)
    # print(a._x_start)
    # print(a._x_end)
    # print(a._y_start)
    # print(a._y_end)
    # print(a._z.shape)
    # print(a._energies)
    # print(a.intensity_as_numpy().shape)

    a.printInfo()
    # a.showEField()


    file_content = numpy.load(filename)

    for k in file_content.keys():
        print(">>>> key: ",k)

    e_field = file_content["e_field"]
    print(e_field.shape)
    coordinates = file_content["coordinates"]
    print(coordinates)
    energies = file_content["energies"]
    print(energies)

    print(coordinates)
    total_intensity = (numpy.abs(e_field[0,:,:,0])**2).sum()
    s = e_field.shape
    x = numpy.linspace(coordinates[0], coordinates[1], s[1])
    y = numpy.linspace(coordinates[2], coordinates[3], s[2])
    dx = 1e0 * (x[1] - x[0])
    dy = 1e0 * (y[1] - y[0])
    print(">>>>",dx,dy)
    plot_image(numpy.abs(e_field[0,:,:,0])**2,1e6*x,1e6*y,show=show,title=title+" Tot Int: %g"%(total_intensity*dx*dy))

    # duplicate test
    # c = NumpyWavefront.fromNumpyArray(e_field, coordinates, energies)
    # print(c)
    # c.printInfo()
    # c.showEField()
if __name__ == "__main__":

    # filename = "/users/srio/Working/paper-hierarchical/CODE-COMSYL/tmp/tmp0_chac_in_mode0.npz"
    # plot_filename(filename=filename,show=True)



    max_mode = 9
    for i in range(max_mode+1):
        filename_in = "/users/srio/Working/paper-hierarchical/CODE-COMSYL/tmp/tmp0_chac_in_mode%d.npz"%i
        filename_out = "/users/srio/Working/paper-hierarchical/CODE-COMSYL/tmp/tmp0_chac_out_mode%d.npz" % i

        plot_filename(filename=filename_in, title="IN %d"%i, show=False)
        if i == max_mode:
            plot_filename(filename=filename_out, title="OUT %d"%i, show=True)
        else:
            plot_filename(filename=filename_out, title="OUT %d"%i, show=False)