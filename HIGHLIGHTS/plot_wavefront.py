

import numpy

from srwlib import *
import sys
from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction
from comsyl.autocorrelation.AutocorrelationFunctionPropagator import AutocorrelationFunctionPropagator
from comsyl.parallel.utils import isMaster, barrier
from comsyl.utils.Logger import log


from comsyl.waveoptics.Wavefront import NumpyWavefront, SRWWavefront

if __name__ == "__main__":


    #
    # wavefronts
    #


    filename = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/propagation/data_free_cs_new_u18_2m_1h_s2.5/cs_new_u18_2m_1h_s2.5_d26.01.wfs.npz"
    filename = "tmp/tmp0_hib3-3302_out.npz"
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
    a.showEField()


    file_content = numpy.load(filename)

    for k in file_content.keys():
        print(">>>> key: ",k)

    e_field = file_content["e_field"]
    print(e_field.shape)
    coordinates = file_content["coordinates"]
    print(coordinates)
    energies = file_content["energies"]
    print(energies)

    c = NumpyWavefront.fromNumpyArray(e_field, coordinates, energies)

