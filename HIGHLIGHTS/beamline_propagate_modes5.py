
import pickle

from comsyl.waveoptics.WOFRYAdapter import ComsylWofryBeamline
from comsyl.waveoptics.SRWAdapter import ComsylSRWBeamline

from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction


def create_beamline_srw_new(load_from_file=None,slit_width=5e-6,slit_height=5e-6,source_offset=0.0,dumpfile=None):

    from srwlib import srwl_uti_proc_is_master, SRWLOptD, SRWLOptL, SRWLOptA, SRWLOptC
    import pickle



    if load_from_file is not None:
        return pickle.load(open(load_from_file,"rb"))

    if not srwl_uti_proc_is_master(): exit()

    ####################################################
    # BEAMLINE

    srw_oe_array = []
    srw_pp_array = []
    drift_before_oe_0 = SRWLOptD(36.0+source_offset)
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


    optBL = SRWLOptC(srw_oe_array, srw_pp_array)


    if dumpfile is not None:
        pickle.dump(ComsylSRWBeamline(optBL), open(dumpfile,"wb"))
        print("File written to disk: ",dumpfile)


    return ComsylSRWBeamline(optBL)

def create_beamline_wofry(load_from_file=None, slit_width=5e-6,slit_height=5e-6,source_offset=0.0,dumpfile=None):

    import pickle

    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen
    from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens
    from wofry.beamline.optical_elements.absorbers.slit import WOSlit
    from syned.beamline.shape import Rectangle
    from wofry.propagator.propagator import PropagationElements, PropagationParameters

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

    # {'shift_half_pixel':1,'magnification_x':1.0,'magnification_y':1.0}

    elements = [beamline_element1,beamline_element2,beamline_element3,beamline_element4]
    handlers = ['FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D']
    specific = [{'shift_half_pixel':1,'magnification_x':8.0,'magnification_y':8.0},
            {'shift_half_pixel':1,'magnification_x':1.0,'magnification_y':1.0},
            {'shift_half_pixel':1,'magnification_x':0.5,'magnification_y':0.35},
            {'shift_half_pixel':1,'magnification_x':0.2,'magnification_y':0.2}]


    propagation_elements = PropagationElements()

    for k in range(len(elements)):
        propagation_elements.add_beamline_element(elements[k], specific[k])

    bl = ComsylWofryBeamline.initialize_from_propagator_elements_object(propagation_elements)

    if dumpfile is not None:
        pickle.dump(bl, open(dumpfile,"wb"))
        print("File written to disk: ",dumpfile)

    return bl


if __name__ == "__main__":


    method = 'WOFRY'

    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++ method =     %s "%(method))
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


    filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS

    autocorrelation_function = AutocorrelationFunction.load(filename)

    #
    # source offfset
    #
    source_position=autocorrelation_function.info().sourcePosition()
    if source_position == "entrance":
        source_offset = autocorrelation_function._undulator.length() * 0.5 #+ 2 * comparer.undulator().periodLength()
        print("Using source position entrance z=%f" % source_offset)
    elif source_position == "center":
        source_offset = 0.0
        print("Using source position center z=%f" % source_offset)
    else:
        raise Exception("Unhandled source position")

    print(">>>>> Using source position center z=%f" % source_offset)

    directory_name = "propagation_EBS_25x25"

    if method == 'SRW':
        #
        # SRW
        #

        comsyl_beamline = create_beamline_srw_new(slit_width=25e-6,slit_height=25e-6,source_offset=0,
                                              dumpfile="bl.p")

        comsyl_beamline.add_undulator_offset(offset=source_offset)
        # print(">>>>>>>>>***************",source_offset,comsyl_beamline.get_native_beamline().arOpt[0].L)

        # comsyl_beamline = create_beamline_srw_new(load_from_file="bl.p")

    elif method == 'WOFRY':
        #
        # WOFRY
        #
        comsyl_beamline = create_beamline_wofry(slit_width=25e-6,slit_height=25e-6,source_offset=0,dumpfile="bl.p")
        print("????OFFSET: ",source_offset)
        print("\n\n????BEFORE: ")
        comsyl_beamline._info()

        comsyl_beamline.add_undulator_offset(source_offset)
        print("\n\n????AFTER: ")
        comsyl_beamline._info()
        # directory_name = "propagation_wofry_EBS_25x25"

    elif method == "FILE":
        comsyl_beamline  = pickle.load(open("bl.p","rb"))

    # print("????",comsyl_beamline._info())
    # print(">>>>>>>>>>>>>>>>>",comsyl_beamline.propagation_code())
    # af_propagated = comsyl_beamline.propagate_af(autocorrelation_function,
    #              directory_name=directory_name,
    #              af_output_file_root="%s/propagated_beamline"%(directory_name),
    #              maximum_mode=2,python_to_be_used="/users/srio/OASYS1.1d/miniconda3/bin/python")
    #
    # print(">>>>>>>>>>>>>>>>>",comsyl_beamline.propagation_code())
