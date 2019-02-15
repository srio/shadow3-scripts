

from srwlib import *

from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction
from comsyl.autocorrelation.AutocorrelationFunctionPropagator import AutocorrelationFunctionPropagator
from comsyl.parallel.utils import isMaster, barrier
from comsyl.utils.Logger import log
# import pickle
from comsyl.waveoptics.WOFRYAdapter import CWBeamline
# import sys


# def create_beamline_srw(distance, undulator, source_offset=0.0):
#
#     # if source_position == "entrance":
#     #     source_offset = undulator.length() * 0.5 #+ 2 * comparer.undulator().periodLength()
#     #     log("Using source position entrance z=%f" % source_offset)
#     # elif source_position == "center":
#     #     source_offset = 0.0
#     #     log("Using source position center z=%f" % source_offset)
#     # else:
#     #     raise Exception("Unhandled source position")
#
#
#     div_x_factor = int(distance) + 1
#     div_y_factor = int(distance) + 1
#
#     optBL = SRWLOptC([SRWLOptD(source_offset+distance)],
#                      [[0, 0, 1.0, 0, 0, div_x_factor, 1,   div_y_factor, 1,    0, 0, 0], [0, 0, 1.0, 0, 0, 1, 0.05/2.0, 1, 0.1, 0, 0, 0]])
#
#
#
#     return optBL

def create_beamline_srw_new(slit_width=5e-6,slit_height=5e-6,source_offset=0.0):

    if not srwl_uti_proc_is_master(): exit()

    ####################################################
    # BEAMLINE

    srw_oe_array = []
    srw_pp_array = []
    drift_before_oe_0 = SRWLOptD(36.0)
    pp_drift_before_oe_0 = [0,0,1.0,0,0,8.0,0.5,10.0,0.5,0,0.0,0.0]

    srw_oe_array.append(drift_before_oe_0)
    srw_pp_array.append(pp_drift_before_oe_0)


    oe_1=SRWLOptL(_Fx=18.0, _Fy=18.0, _x=0.0, _y=0.0)
    pp_oe_1 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_1)
    srw_pp_array.append(pp_oe_1)

    drift_before_oe_2 = SRWLOptD(32.0)
    pp_drift_before_oe_2 = [0,0,1.0,0,0,0.4,1.0,0.25,1.0,0,0.0,0.0]

    srw_oe_array.append(drift_before_oe_2)
    srw_pp_array.append(pp_drift_before_oe_2)

    oe_2=SRWLOptA(_shape='r',
                   _ap_or_ob='a',
                   _Dx=slit_width,
                   _Dy=slit_height,
                   _x=0.0,
                   _y=0.0)
    pp_oe_2 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

    srw_oe_array.append(oe_2)
    srw_pp_array.append(pp_oe_2)

    drift_before_oe_3 = SRWLOptD(2.0)
    pp_drift_before_oe_3 = [0,0,1.0,0,0,0.3,5.0,0.25,2.0,0,0.0,0.0]

    srw_oe_array.append(drift_before_oe_3)
    srw_pp_array.append(pp_drift_before_oe_3)



    ####################################################
    # PROPAGATION

    optBL = SRWLOptC(srw_oe_array, srw_pp_array)




    return optBL

def create_beamline_wofry(load_from_file=None, slit_width=5e-6,slit_height=5e-6,source_offset=0.0,dumpfile=None):

    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen
    from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens
    from wofry.beamline.optical_elements.absorbers.slit import WOSlit
    from syned.beamline.shape import Rectangle
    from wofry.propagator.propagator import PropagationElements, PropagationParameters
    from syned.beamline.beamline import Beamline

    from syned.storage_ring.empty_light_source import EmptyLightSource

    # from comsyl.waveoptics.WofrySuperBeamline import WofrySuperBeamline
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


    propagation_elements = PropagationElements()

    for k in range(len(elements)):
        propagation_elements.add_beamline_element(elements[k], specific[k])

    bl = CWBeamline.initialize_from_propagator_elements_object(propagation_elements)

    if dumpfile is not None:
        pickle.dump(bl, open(dumpfile,"wb"))
        print("File written to disk: ",dumpfile)

    return bl



def propagate_modes_srw(beamline, directory_name,maximum_mode=None):

    propagator = AutocorrelationFunctionPropagator(beamline)

    if maximum_mode is None:
        mode_distribution=autocorrelation_function.modeDistribution()
        maximum_mode = mode_distribution[abs(mode_distribution)>0.00005].shape[0]

    propagator.setMaximumMode(maximum_mode)
    data_directory = "%s/beamline_%s" % (directory_name, af_name)

    if isMaster():
        if not os.path.exists(directory_name):
            os.mkdir(directory_name)
        if not os.path.exists(data_directory):
            os.mkdir(data_directory)
    barrier()

    propagated_filename = "%s/%s.npz" % (data_directory, af_name)
    af = propagator.propagate(autocorrelation_function, propagated_filename,
                              python_to_be_used="/users/srio/OASYS1.1d/miniconda3/bin/python")
    af.save("%s/beamline_%s.npz" % (directory_name, af_name))


if __name__ == "__main__":

    import sys

    if isMaster():
        print("This is the name of the script: ", sys.argv[0] )
        print("Number of arguments: ", len(sys.argv) )
        print("The arguments are: " , str(sys.argv) )

        if len(sys.argv) == 1:

            method = 'WOFRY'
            ring = "EBS" # "LB"  "HB"
            VERTICAL = "25"
            HORIZONTAL = "25"
        else:
            method =     sys.argv[1]
            ring =       sys.argv[2]
            VERTICAL =   sys.argv[3]
            HORIZONTAL = sys.argv[4]


        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("++++ method =     %s "%(method))
        print("++++ ring =       %s "%(ring))
        print("++++ VERTICAL =   %s "%(VERTICAL))
        print("++++ HORIZONTAL = %s "%(HORIZONTAL))
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    if ring == "EBS":
        filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    elif ring == "HB":
        filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s2.0.npz"
    elif ring == "LB":
        filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npz" # OK LB


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


    slit_height = float(VERTICAL)*1e-6
    slit_width = float(HORIZONTAL)*1e-6

    if method == 'SRW':
        #
        # SRW
        #
        af_name = filename.split("/")[-1].replace(".npz", "")
        directory_name = "propagation_srw_%s_%sx%s"%(ring,VERTICAL,HORIZONTAL)
        distance = 26.0 # float(sys.argv[1])
        undulator = autocorrelation_function._undulator
        # beamlineSRW = create_beamline_srw(distance, undulator, source_offset=source_offset)
        beamlineSRW = create_beamline_srw_new(slit_width=slit_width,slit_height=slit_height,source_offset=source_offset)
        propagate_modes_srw(beamlineSRW, directory_name, maximum_mode=2)
    elif method == 'WOFRY':
        #
        # WOFRY
        #
        # beamlineWOFRY = create_beamline_wofry(slit_width=slit_width,slit_height=slit_height,source_offset=source_offset,dumpfile="tmp.p")
        beamlineWOFRY = create_beamline_wofry(load_from_file="tmp.p")

        print(beamlineWOFRY)
        directory_name = "propagation_wofry_%s_%sx%s"%(ring,VERTICAL,HORIZONTAL)
        af_propagated = beamlineWOFRY.propagate_af(autocorrelation_function,
                     directory_name=directory_name,
                     af_output_file_root="%s/propagated_beamline"%(directory_name),
                     maximum_mode=2,python_to_be_used="/users/srio/OASYS1.1d/miniconda3/bin/python")
