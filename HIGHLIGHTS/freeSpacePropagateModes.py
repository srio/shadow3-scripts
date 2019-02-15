__author__ = 'mglass'

from srwlib import *
import sys
from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction
from comsyl.autocorrelation.AutocorrelationFunctionPropagator import AutocorrelationFunctionPropagator
from comsyl.parallel.utils import isMaster, barrier
from comsyl.utils.Logger import log


def createBeamlinePS(distance, undulator, source_position):

    if source_position == "entrance":
        source_offset = undulator.length() * 0.5 #+ 2 * comparer.undulator().periodLength()
        log("Using source position entrance z=%f" % source_offset)
    elif source_position == "center":
        source_offset = 0.0
        log("Using source position center z=%f" % source_offset)
    else:
        raise Exception("Unhandled source position")


    div_x_factor = int(distance) + 1
    div_y_factor = int(distance) + 1

    optBL = SRWLOptC([SRWLOptD(source_offset+distance)],
                     [[0, 0, 1.0, 0, 0, div_x_factor, 1,   div_y_factor, 1,    0, 0, 0], [0, 0, 1.0, 0, 0, 1, 0.05/2.0, 1, 0.1, 0, 0, 0]])



    return optBL


def propagateModes(distance, filename, directory_name,maximum_mode=None):

    af_name = filename.split("/")[-1].replace(".npz", "")

    autocorrelation_function = AutocorrelationFunction.load(filename)

    undulator = autocorrelation_function._undulator

    beamline = createBeamlinePS(distance, undulator, source_position=autocorrelation_function.info().sourcePosition())

    propagator = AutocorrelationFunctionPropagator(beamline)

    if maximum_mode is None:
        mode_distribution=autocorrelation_function.modeDistribution()
        maximum_mode = mode_distribution[abs(mode_distribution)>0.00005].shape[0]

    propagator.setMaximumMode(maximum_mode)
    data_directory = "%s/data_free_%s" % (directory_name, af_name)

    if isMaster():
        if not os.path.exists(data_directory):
            os.mkdir(data_directory)
    barrier()


    propagated_filename = "%s/%s_d%.1f.npz" % (data_directory, af_name, distance)
    af = propagator.propagate(autocorrelation_function, propagated_filename)
    af.save("%s/free_prop_%s_d%.1f.npz" % (directory_name, af_name, distance))

if __name__ == "__main__":
    # if len(sys.argv) <= 2:
    #     print("Need distance and filename")
    #     exit()

    filename_ebs = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    # filename_lb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy" # OK LB
    # filename_hb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s2.0.npy"


    distance = 26.0 # float(sys.argv[1])
    filename = filename_ebs # sys.argv[2]
    directory_name = "propagation"

    propagateModes(distance, filename, directory_name, maximum_mode=50)
