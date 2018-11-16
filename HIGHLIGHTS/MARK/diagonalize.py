

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
    # AF
    #

    cases = ["15x15","10x10","5x5"]

    for case in cases:
        filename = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/NICE/propagation_wofry_%s/beamline_cs_new_u18_2m_1h_s2.5.npz"%case

        af_name = filename.split("/")[-1].replace(".npz", "")

        autocorrelation_function = AutocorrelationFunction.load(filename)

        print(autocorrelation_function.eigenvalues())

        autocorrelation_function.diagonalizeModes(10)

        print(autocorrelation_function.eigenvalues())

        autocorrelation_function.save("rediagonalized_%s"%case)


