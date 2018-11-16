#
# runs using functions to propagate individual elements" ans use wofry wavefronts
#


from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
# from CompactAFReader import CompactAFReader

import numpy


from srxraylib.plot.gol import plot_image, plot


# from plot_color import plot_with_transparency_one

import pylab as plt
from matplotlib.colors import Normalize, ListedColormap
import matplotlib.patches as patches



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


def propagator_initialize():

    from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
    from wofry.propagator.propagator import PropagationManager

    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoomXY2D())
    except:
        pass

    return propagator

def propagate_screen1(input_wavefront,magnification_x=8.0,magnification_y=8.0):


    from wofry.propagator.propagator import PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen

    optical_element = WOScreen()
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()

    beamline_element = BeamlineElement(optical_element=optical_element,coordinates=ElementCoordinates(
        p=36.0,q=0,angle_radial=numpy.radians(0.0),angle_azimuthal=numpy.radians(0.0)))

    propagation_elements.add_beamline_element(beamline_element)

    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
    propagation_parameters.set_additional_parameters('magnification_x', magnification_x)
    propagation_parameters.set_additional_parameters('magnification_y', magnification_y)

    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,handler_name='FRESNEL_ZOOM_XY_2D')
    return output_wavefront


def propagate_lens(input_wavefront,magnification_x=5.0,magnification_y=0.1):

    from wofry.propagator.propagator import  PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens


    optical_element = WOIdealLens(name='',focal_x=18.000000,focal_y=18.000000)

    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=0.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),    propagation_elements = propagation_elements)
    #
    propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
    propagation_parameters.set_additional_parameters('magnification_x', magnification_x)
    propagation_parameters.set_additional_parameters('magnification_y', magnification_y)
    #
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_XY_2D')

    return output_wavefront

def propagate_slit(input_wavefront,magnification_x=0.50,magnification_y=0.35,slit_width=5e-6,slit_height=5e-6):

    from wofry.propagator.propagator import PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from syned.beamline.shape import Rectangle

        #
    # define current oe
    #
    boundary_shape = Rectangle(x_left=-0.5*slit_width,x_right=0.5*slit_width,y_bottom=-0.5*slit_height,y_top=0.5*slit_height)

    from wofry.beamline.optical_elements.absorbers.slit import WOSlit
    optical_element = WOSlit(boundary_shape=boundary_shape)

    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=32.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
    propagation_parameters.set_additional_parameters('magnification_x', magnification_x)
    propagation_parameters.set_additional_parameters('magnification_y', magnification_y)
    #
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_XY_2D')

    return output_wavefront


def propagate_screen2(input_wavefront,magnification_x=0.2,magnification_y=0.2):


    from wofry.propagator.propagator import PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen

    optical_element = WOScreen()

    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
            coordinates=ElementCoordinates(p=2.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),    propagation_elements = propagation_elements)
    #self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
    propagation_parameters.set_additional_parameters('magnification_x', magnification_x)
    propagation_parameters.set_additional_parameters('magnification_y', magnification_y)

    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_XY_2D')
    return output_wavefront


if __name__ == "__main__":


    filename_ebs = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    filename_lb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy" # OK LB
    filename_hb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s2.0.npy"

    #
    # load CSD
    #

    af  = CompactAFReader.initialize_from_file(filename_ebs)
    x = af.x_coordinates()
    y = af.y_coordinates()
    cumulated_occupation = af.cumulated_occupation_array()
    occupation = af.occupation_array()

    # plot(numpy.arange(cumulated_occupation.size),cumulated_occupation)


    # print(af.info())
    propagator = propagator_initialize()

    intensity_full = af.total_intensity()
    print(">>>>>>>>>>>>> intensity (full)",intensity_full, " from modes: ",af.total_intensity_from_modes())


    f = open("propagation.txt",'w')
    f.close()

    slitV = [30,25,20,15,10,5, 5]
    slitH = [25,25,20,15,10,10,5]

    for j in range(2): #len(slitV)):
        intensity_accumulated = 0.0
        intensity_accumulated2 = 0.0
        for i in range(5): #af.number_modes()):


            # eigen_m = af.mode(i)
            # eigen_v = af.eigenvalue(i).real
            # wf = create_wofry_wavefront(af,i)
            # deltas = (wf.get_coordinate_x()[1]-wf.get_coordinate_x()[0])*(wf.get_coordinate_y()[1]-wf.get_coordinate_y()[0])
            # print(i,eigen_v,(eigen_m.conjugate()*eigen_m).sum()*deltas,"<><> %g"%(wf.get_intensity().sum()*deltas))




            wf0 = create_wofry_wavefront(af,i)

            wf1 = propagate_screen1(wf0)

            wf2 = propagate_lens(wf1)

            wf3 = propagate_slit(wf2,magnification_x=0.50,magnification_y=0.35,
                                 # slit_width=1,slit_height=1)
                                 slit_width=slitH[j]*1e-6,slit_height=slitV[j]*1e-6)

            wf4 = propagate_screen2(wf3)


            intensity_accumulated += wf4.get_integrated_intensity()

            # wf_list = [wf0,wf1,wf2,wf3,wf4]
            # for wf in wf_list:
            #     deltas = (wf.get_coordinate_x()[1]-wf.get_coordinate_x()[0])*(wf.get_coordinate_y()[1]-wf.get_coordinate_y()[0])
            #     print(">>>>>>>>>>>>>> mode: %d intensity: %g"%(i,wf.get_intensity().sum()*deltas))

            # plot_image(wf4.get_intensity(),1e6*wf4.get_coordinate_x(),1e6*wf4.get_coordinate_y(),title="mode %d propagated to final oe"%i)

            ratio = intensity_accumulated/intensity_full
            print("up to mode: %d intensity(full): %g intensity(cut): %g ratio: %g "%(i,intensity_full,intensity_accumulated,ratio,))


        f = open("propagation.txt",'a')
        f.write('"%d x %d" %f \n'%(slitV[j],slitH[j],ratio))
        f.close