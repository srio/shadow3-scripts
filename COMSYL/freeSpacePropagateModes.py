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
    file_out = "%s/free_prop_%s_d%.1f.npz" % (directory_name, af_name, distance)
    af.save(file_out)
    print("File written to disk: "+file_out)

if __name__ == "__main__":
    # if len(sys.argv) <= 2:
    #     print("Need distance and filename")
    #     exit()
    # distance = float(sys.argv[1])
    # filename = sys.argv[2]
    directory_name = "propagation"

    filename = "/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npz"

    distance = 10.0

    propagateModes(distance, filename, directory_name, maximum_mode=3)