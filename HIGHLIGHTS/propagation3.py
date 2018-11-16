#
# runs using "WofrySuperBeamline" adapted for numpy wavefronts ans pickle
#



import numpy


# from srxraylib.plot.gol import plot_image, plot
# from wofry.propagator.propagator import PropagationElements, PropagationParameters
# from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D







def create_wofry_wavefront(af,index,weight_with_eigenvalue=True):

    from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D

    z_array = af.mode(index)

    if weight_with_eigenvalue:
        z_array *= numpy.sqrt( (af.eigenvalue(index).real ))

    wavefront = GenericWavefront2D.initialize_wavefront_from_arrays(x_array=af.x_coordinates(),
                                                                          y_array=af.y_coordinates(),
                                                                          z_array=z_array,
                                                                          )

    wavefront.set_photon_energy(af.photon_energy())

    return wavefront


def create_wofry_beamline(load_from_file=None, slit_width=5e-6,slit_height=5e-6):

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
            coordinates=ElementCoordinates(p=36.0,q=0))

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


if __name__ == "__main__":

    from comsyl.waveoptics.WofrySuperBeamline import WofrySuperBeamline
    from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
    # from CompactAFReader import CompactAFReader
    from wofry.propagator.propagator import PropagationManager
    from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D


    # initialize propagator
    mypropagator = PropagationManager.Instance()
    try:
        mypropagator.add_propagator(FresnelZoomXY2D())
    except:
        print("May be you alreay initialized propagator and stored FresnelZoomXY2D")


    filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    # filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy" # OK LB
    # filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s2.0.npy"

    #
    # load CSD
    #

    af  = CompactAFReader.initialize_from_file(filename)
    x = af.x_coordinates()
    y = af.y_coordinates()
    cumulated_occupation = af.cumulated_occupation_array()
    occupation = af.occupation_array()

    intensity_full = af.total_intensity()

    # plot(numpy.arange(cumulated_occupation.size),cumulated_occupation)


    # customize beamline
    slit_width=5e-6
    slit_height=5e-6

    slitV = [30,25,20,15,10,5, 5]
    slitH = [25,25,20,15,10,10,5]



    #
    # loop
    #
    f = open("propagation.txt",'w')
    f.close()


    for j in range(1): #len(slitV)):

        slit_width = slitH[j]*1e-6
        slit_height = slitV[j]*1e-6

        intensity_accumulated = 0.0
        intensity_accumulated2 = 0.0

        for i in range(3): #af.number_modes()):


            BEAMLINE = create_wofry_beamline(load_from_file=None, #"./BEAMLINE.p",
                slit_width=slitH[j]*1e-6,slit_height=slitV[j]*1e-6)


            wofry_wavefront = create_wofry_wavefront(af,i,weight_with_eigenvalue=True)
            wf_list = BEAMLINE.propagate(wofry_wavefront,mypropagator)
            # wf_list = WofrySuperBeamline.propagate_numpy_wavefront(
            #     # "tmp/tmp0_hib3-3302_out.npz",
            #     "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/tmp-working/tmp0_hib3-3302_IN.npz",
            #     "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/tmp-working/tmp0_hib3-3302_OUT.npz",
            #     BEAMLINE,mypropagator)


            intensity_accumulated += wf_list[-1].get_integrated_intensity()

            ratio = intensity_accumulated/intensity_full
            print("slit H:%d, V:%d, up to mode: %d intensity(full): %g intensity(cut): %g ratio: %g "%(
                slitH[j],slitV[j],i,intensity_full,intensity_accumulated,ratio,))


            print("\n\n\n\n\n\n")
        f = open("propagation.txt",'a')
        f.write('"%d x %d" %f \n'%(slitV[j],slitH[j],ratio))
        f.close


