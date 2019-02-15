__author__ = 'mglass'

from srwlib import *

from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction
from comsyl.autocorrelation.AutocorrelationFunctionPropagator import AutocorrelationFunctionPropagator
from comsyl.parallel.utils import isMaster, barrier
from comsyl.utils.Logger import log
# import pickle
# from comsyl.waveoptics.WofrySuperBeamline import WofrySuperBeamline
# import sys


def create_beamline_srw(distance, undulator, source_offset=0.0):

    # if source_position == "entrance":
    #     source_offset = undulator.length() * 0.5 #+ 2 * comparer.undulator().periodLength()
    #     log("Using source position entrance z=%f" % source_offset)
    # elif source_position == "center":
    #     source_offset = 0.0
    #     log("Using source position center z=%f" % source_offset)
    # else:
    #     raise Exception("Unhandled source position")


    div_x_factor = int(distance) + 1
    div_y_factor = int(distance) + 1

    optBL = SRWLOptC([SRWLOptD(source_offset+distance)],
                     [[0, 0, 1.0, 0, 0, div_x_factor, 1,   div_y_factor, 1,    0, 0, 0], [0, 0, 1.0, 0, 0, 1, 0.05/2.0, 1, 0.1, 0, 0, 0]])



    return optBL

# def create_beamline_wofry(load_from_file=None):
#     if load_from_file is not None:
#         return pickle.load(open(load_from_file,"rb"))
def create_beamline_wofry(load_from_file=None, slit_width=5e-6,slit_height=5e-6,source_offset=0.0):

    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen
    from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens
    from wofry.beamline.optical_elements.absorbers.slit import WOSlit
    from syned.beamline.shape import Rectangle

    from comsyl.waveoptics.WofrySuperBeamline import WofrySuperBeamline
    import pickle
    #
    # get the list of beamline elements and the propagator specific parameters
    #
    if load_from_file is not None:
        return pickle.load(open(load_from_file,"rb"))


    # first screen
    beamline_element1 = BeamlineElement(
            optical_element=WOScreen(),
            coordinates=ElementCoordinates(p=36.0+source_offset,q=0))

    # lens
    beamline_element2 = BeamlineElement(
        optical_element=WOIdealLens(name='',focal_x=18.000000,focal_y=18.000000),
        coordinates=ElementCoordinates(p=0,q=0))

    # slit


    boundary_shape = Rectangle(x_left=-0.5*slit_width,x_right=0.5*slit_width,y_bottom=-0.5*slit_height,y_top=0.5*slit_height)
    beamline_element3 = BeamlineElement(
        optical_element=WOSlit(boundary_shape=boundary_shape),
        coordinates=ElementCoordinates(p=32,q=0))

    # screen
    beamline_element4 = BeamlineElement(
        optical_element=WOScreen(),
        coordinates=ElementCoordinates(p=2,q=0))

    {'shift_half_pixel':1,'magnification_x':1.0,'magnification_y':1.0}

    elements = [beamline_element1,beamline_element2,beamline_element3,beamline_element4]
    handlers = ['FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D']
    specific = [{'shift_half_pixel':1,'magnification_x':8.0,'magnification_y':8.0},
            {'shift_half_pixel':1,'magnification_x':1.0,'magnification_y':1.0},
            {'shift_half_pixel':1,'magnification_x':0.5,'magnification_y':0.35},
            {'shift_half_pixel':1,'magnification_x':0.2,'magnification_y':0.2}]


    # create beamline
    BEAMLINE = WofrySuperBeamline()

    for k in range(len(elements)):
        BEAMLINE.add_element(beamline_element=elements[k],
                             propagator_handler=handlers[k],
                             propagator_specific_parameters=specific[k])

    # pickle.dump(BEAMLINE, open("./BEAMLINE.p","wb"))


    return BEAMLINE


def propagate_modes_srw(beamline, directory_name,maximum_mode=None):

    propagator = AutocorrelationFunctionPropagator(beamline)

    if maximum_mode is None:
        mode_distribution=autocorrelation_function.modeDistribution()
        maximum_mode = mode_distribution[abs(mode_distribution)>0.00005].shape[0]

    propagator.setMaximumMode(maximum_mode)
    data_directory = "%s/beamline_%s" % (directory_name, af_name)

    if isMaster():
        if not os.path.exists(data_directory):
            os.mkdir(data_directory)
    barrier()

    #propagated_filename = "%s/%s_d%.1f.npz" % (data_directory, af_name, distance)

    propagated_filename = "%s/%s.npz" % (data_directory, af_name)
    af = propagator.propagate(autocorrelation_function, propagated_filename,
                              python_to_be_used="/users/srio/OASYS1.1/miniconda3/bin/python")
    af.save("%s/beamline_%s.npz" % (directory_name, af_name))

def propagate_modes_wofry(beamline, directory_name,maximum_mode=None):

    propagator = AutocorrelationFunctionPropagator(beamline)

    if maximum_mode is None:
        mode_distribution=autocorrelation_function.modeDistribution()
        maximum_mode = mode_distribution[abs(mode_distribution)>0.00005].shape[0]

    propagator.setMaximumMode(maximum_mode)
    data_directory = "%s/beamline_%s" % (directory_name, af_name)

    if isMaster():
        if not os.path.exists(data_directory):
            os.mkdir(data_directory)
    barrier()


    propagated_filename = "%s/%s.npz" % (data_directory, af_name)
    af = propagator.propagate(autocorrelation_function, propagated_filename,method="WOFRY",
                              python_to_be_used="/users/srio/OASYS1.1/miniconda3/bin/python")
    af.save("%s/beamline_%s.npz" % (directory_name, af_name))


if __name__ == "__main__":
    import os

    os.system('ls')
    # Load initial modes
    filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    af_name = filename.split("/")[-1].replace(".npz", "")
    autocorrelation_function = AutocorrelationFunction.load(filename)

    #
    # source offfset
    #
    source_position=autocorrelation_function.info().sourcePosition()
    if source_position == "entrance":
        source_offset = autocorrelation_function._undulator.length() * 0.5 #+ 2 * comparer.undulator().periodLength()
        log("Using source position entrance z=%f" % source_offset)
    elif source_position == "center":
        source_offset = 0.0
        log("Using source position center z=%f" % source_offset)
    else:
        raise Exception("Unhandled source position")

    print(">>>>> Using source position center z=%f" % source_offset)




    method = 'WOFRY'

    if method == 'SRW':
        #
        # SRW
        #
        directory_name = "propagation_srw"
        distance = 26.0 # float(sys.argv[1])
        undulator = autocorrelation_function._undulator
        beamlineSRW = create_beamline_srw(distance, undulator, source_offset=source_offset)
        propagate_modes_srw(beamlineSRW, directory_name, maximum_mode=1)
    elif method == 'WOFRY':
        #
        # WOFRY
        #
        directory_name = "propagation_wofry"
        beamlineWOFRY = create_beamline_wofry(slit_width=25e-6,slit_height=25e-6,source_offset=source_offset) #load_from_file="/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/BEAMLINE.p")
        propagate_modes_wofry(beamlineWOFRY, directory_name,maximum_mode=2)


        # os.system()



