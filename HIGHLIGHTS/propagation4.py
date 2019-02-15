#
# runs using modified wofry propagator and adapted for numpy wavefronts ans pickle
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

def create_beamline_wofry(load_from_file=None, slit_width=5e-6,slit_height=5e-6,source_offset=0.0):

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


    # create beamline
    # BEAMLINE = WofrySuperBeamline()

    propagation_elements = PropagationElements()

    for k in range(len(elements)):
        propagation_elements.add_beamline_element(elements[k], specific[k])

    # pickle.dump(BEAMLINE, open("./BEAMLINE.p","wb"))

    # return propagation_elements


    #
    # beamline
    #

    # beamline1 = Beamline(beamline_elements_list=[])
    beamline1 = Beamline()
    print(">>>>????",beamline1._beamline_elements_list)
    beamline1.set_light_source(EmptyLightSource())

    print(">>> beamline before",beamline1.get_beamline_elements_number())
    for k in range(len(elements)):
        beamline1.append_beamline_element(elements[k]) # , specific[k])
    print(">>> beamline after",beamline1.get_beamline_elements_number())

    # beamline1.append_beamline_element(BeamlineElement(screen1,coordinates=ElementCoordinates(100.0)))

    return beamline1, specific


if __name__ == "__main__":

    from comsyl.waveoptics.WofrySuperBeamline import WofrySuperBeamline
    from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
    # from CompactAFReader import CompactAFReader
    from wofry.propagator.propagator import PropagationManager
    from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D

    from wofry.propagator.propagator import PropagationElements, PropagationParameters


    # initialize propagator
    mypropagator = PropagationManager.Instance()
    try:
        mypropagator.add_propagator(FresnelZoomXY2D())
    except:
        print("May be you alreay initialized propagator and stored FresnelZoomXY2D")



    #
    # beamline, specific = create_beamline_wofry()
    #
    # propagation_elements = PropagationElements()
    # for i,element in enumerate(beamline.get_beamline_elements()):
    #     propagation_elements.add_beamline_element(beamline_element=element,
    #                                               element_parameters=specific[i])





    # wofry_wavefront = create_wofry_wavefront()
    # print(wofry_wavefront)
    # propagation_parameters = PropagationParameters(wavefront=wofry_wavefront.duplicate(),
    #                                                propagation_elements = propagation_elements)


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




            wofry_wavefront = create_wofry_wavefront(af,i,weight_with_eigenvalue=True)


            # propagation_elements = create_beamline_wofry(load_from_file=None, #"./BEAMLINE.p",
            #     slit_width=slitH[j]*1e-6,slit_height=slitV[j]*1e-6)
            # propagation_parameters = PropagationParameters(wavefront=wofry_wavefront.duplicate(),propagation_elements = propagation_elements)


            beamline, specific = create_beamline_wofry(load_from_file=None, #"./BEAMLINE.p",
                slit_width=slitH[j]*1e-6,slit_height=slitV[j]*1e-6)


            propagation_elements = PropagationElements()
            for i in range(beamline.get_beamline_elements_number()):
                print("ADDING ELEMENT: ",i,beamline.get_beamline_elements_number(),len(specific))
                propagation_elements.add_beamline_element(beamline_element=beamline.get_beamline_element_at(i),
                                                         element_parameters=specific[i])

            propagation_parameters = PropagationParameters(wavefront=wofry_wavefront.duplicate(),propagation_elements = propagation_elements)

            wf = mypropagator.do_propagation(propagation_parameters=propagation_parameters,handler_name='FRESNEL_ZOOM_XY_2D')


            # wf_list = WofrySuperBeamline.propagate_numpy_wavefront(
            #     # "tmp/tmp0_hib3-3302_out.npz",
            #     "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/tmp-working/tmp0_hib3-3302_IN.npz",
            #     "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/tmp-working/tmp0_hib3-3302_OUT.npz",
            #     BEAMLINE,mypropagator)


            intensity_accumulated += wf.get_integrated_intensity()

            ratio = intensity_accumulated/intensity_full
            print("slit H:%d, V:%d, up to mode: %d intensity(full): %g intensity(cut): %g ratio: %g "%(
                slitH[j],slitV[j],i,intensity_full,intensity_accumulated,ratio,))


            print("\n\n\n\n\n\n")
        f = open("propagation.txt",'a')
        f.write('"%d x %d" %f \n'%(slitV[j],slitH[j],ratio))
        f.close
        print("File written to disk: propagation.txt")


