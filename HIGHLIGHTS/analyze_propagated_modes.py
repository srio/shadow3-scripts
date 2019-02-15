__author__ = 'mglass'

# # from srwlib import *
# import sys
# from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction
from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader


# from comsyl.autocorrelation.AutocorrelationFunctionPropagator import AutocorrelationFunctionPropagator
# from comsyl.parallel.utils import isMaster, barrier
# from comsyl.utils.Logger import log
# import pickle
# from comsyl.waveoptics.WofrySuperBeamline import WofrySuperBeamline


import numpy


if __name__ == "__main__":

    # Load initial modes
    filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    af_name = filename.split("/")[-1].replace(".npz", "")
    autocorrelation_function1 = CompactAFReader.initialize_from_file(filename)
    intensity_full1 = autocorrelation_function1.total_intensity()


    # Load final modes
    filename = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/propagation_wofry/beamline_cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    filename = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/NICE/propagation_wofry_25x25/beamline_cs_new_u18_2m_1h_s2.5.npz"
    af_name = filename.split("/")[-1].replace(".npz", "")
    autocorrelation_function2 = CompactAFReader.initialize_from_file(filename)
    intensity_full2 = autocorrelation_function2.total_intensity()

    eigenvalues1 = autocorrelation_function1.eigenvalues()
    eigenvalues2 = autocorrelation_function2.eigenvalues()

    occupation1 = autocorrelation_function1.occupation_array()
    occupation2 = autocorrelation_function2.occupation_array()

    cumulated_occupation1 = autocorrelation_function1.cumulated_occupation_array()
    cumulated_occupation2 = autocorrelation_function2.cumulated_occupation_array()

    print("Full intensity: before: %g, after: %g, ratio: %g"%(intensity_full1,intensity_full2,intensity_full2/intensity_full1))

    tot1 = autocorrelation_function1.total_intensity()
    tot2 = autocorrelation_function2.total_intensity()

    for i in range(15): # eigenvalues2.size):
        # print(i,eigenvalues1[i],eigenvalues2[i],numpy.abs(eigenvalues2[i]/eigenvalues1[i]))
        # print(i,cumulated_occupation1[i],cumulated_occupation2[i],numpy.abs(cumulated_occupation2[i]/cumulated_occupation1[i]))

        # a1 = numpy.abs(eigenvalues1[i]/eigenvalues1.sum())
        # a2 = numpy.abs(eigenvalues2[i]/eigenvalues2.sum())

        # a1 = numpy.abs(eigenvalues1[i]/tot1)
        # a2 = numpy.abs(eigenvalues2[i]/tot2)
        # print(i,a1,a2,a2/a1)

        a1 = occupation1[i]
        a2 = occupation2[i]

        m1 = autocorrelation_function1.mode(i)
        m2 = autocorrelation_function2.mode(i)

        print(">>>>",autocorrelation_function1.mode_intensity(i),autocorrelation_function2.mode_intensity(i))
        # print(">>",a1,a2,a2/a1,)

    print(">>>> size",eigenvalues1.size,eigenvalues2.size)
    print(">>>> intensities",tot1,tot2)

    #
    # rediagonalize
    #

    # print("BEFORE: ",autocorrelation_function2.eigenvalues())
    # # print("BEFORE RATIO IND: ",autocorrelation_function2.eigenvalues()[0:14]/autocorrelation_function1.eigenvalues()[0:14])
    # print("BEFORE (NORM): ",autocorrelation_function1.eigenvalues()/intensity_full1)
    # print("BEFORE (NORM): ",autocorrelation_function2.eigenvalues()/intensity_full1)

    # autocorrelation_function2._af.diagonalizeModes(10)
    # print("AFTER: ",autocorrelation_function2.eigenvalues())









