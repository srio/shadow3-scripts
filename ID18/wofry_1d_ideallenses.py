#
# Import section
#
import numpy

from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

from wofryimpl.propagator.propagators1D.fresnel import Fresnel1D
from wofryimpl.propagator.propagators1D.fresnel_convolution import FresnelConvolution1D
from wofryimpl.propagator.propagators1D.fraunhofer import Fraunhofer1D
from wofryimpl.propagator.propagators1D.integral import Integral1D
from wofryimpl.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
from wofryimpl.propagator.propagators1D.fresnel_zoom_scaling_theorem import FresnelZoomScaling1D

from srxraylib.plot.gol import plot, plot_image


from syned.beamline.shape import *



#
#
#


#
# ===== Example of python code to create propagate current element =====
#

#
# Import section
#
# import numpy
# from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
# from syned.beamline.beamline_element import BeamlineElement
# from syned.beamline.element_coordinates import ElementCoordinates
# from wofry.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
# from srxraylib.plot.gol import set_qt


def get_wavefront_intensity_fwhm(wf):
    from oasys.util.oasys_util import get_fwhm
    fwhm, quote, coordinates = get_fwhm(wf.get_intensity(), wf.get_abscissas())
    return fwhm

def get_wavefront_intensity_I0(wf):
    I = wf.get_intensity()
    return I[I.size // 2]


# def first_lens_focus(q=50.0):
#     #
#     # Import section
#     #
#     # import numpy
#     #
#     # from syned.beamline.beamline_element import BeamlineElement
#     # from syned.beamline.element_coordinates import ElementCoordinates
#     # from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
#     #
#     # from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
#     #
#     # from wofryimpl.propagator.propagators1D.fresnel import Fresnel1D
#     # from wofryimpl.propagator.propagators1D.fresnel_convolution import FresnelConvolution1D
#     # from wofryimpl.propagator.propagators1D.fraunhofer import Fraunhofer1D
#     # from wofryimpl.propagator.propagators1D.integral import Integral1D
#     # from wofryimpl.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
#     # from wofryimpl.propagator.propagators1D.fresnel_zoom_scaling_theorem import FresnelZoomScaling1D
#     #
#     # from srxraylib.plot.gol import plot, plot_image
#
#     ##########  SOURCE ##########
#
#     #
#     # create output_wavefront
#     #
#     #
#     output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
#                                                                           number_of_points=5000)
#     output_wavefront.set_photon_energy(10000)
#     output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=0, shift=0, beta=0.0922395)
#
#     # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')
#
#     ##########  OPTICAL SYSTEM ##########
#
#     ##########  OPTICAL ELEMENT NUMBER 1 ##########
#
#     input_wavefront = output_wavefront.duplicate()
#     from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
#
#     optical_element = WOScreen1D()
#
#     # drift_before 35 m
#     #
#     # propagating
#     #
#     #
#     propagation_elements = PropagationElements()
#     beamline_element = BeamlineElement(optical_element=optical_element,
#                                        coordinates=ElementCoordinates(p=35.000000, q=0.000000,
#                                                                       angle_radial=numpy.radians(0.000000),
#                                                                       angle_azimuthal=numpy.radians(0.000000)))
#     propagation_elements.add_beamline_element(beamline_element)
#     propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
#     # self.set_additional_parameters(propagation_parameters)
#     #
#     propagation_parameters.set_additional_parameters('magnification_x', 8.0)
#     #
#     propagator = PropagationManager.Instance()
#     try:
#         propagator.add_propagator(FresnelZoom1D())
#     except:
#         pass
#     output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
#                                                  handler_name='FRESNEL_ZOOM_1D')
#
#     #
#     # ---- plots -----
#     #
#     # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 1')
#
#     ##########  OPTICAL ELEMENT NUMBER 2 ##########
#
#     input_wavefront = output_wavefront.duplicate()
#     # from syned.beamline.shape import *
#     boundary_shape = Rectangle(-1.25e-05, 1.25e-05, -1.25e-05, 1.25e-05)
#     from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
#     optical_element = WOSlit1D(boundary_shape=boundary_shape)
#
#     # no drift in this element
#     output_wavefront = optical_element.applyOpticalElement(input_wavefront)
#
#     #
#     # ---- plots -----
#     #
#     # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 2')
#
#     ##########  OPTICAL ELEMENT NUMBER 3 ##########
#
#     input_wavefront = output_wavefront.duplicate()
#     from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
#
#     optical_element = WOScreen1D()
#
#     # drift_before 20 m
#     #
#     # propagating
#     #
#     #
#     propagation_elements = PropagationElements()
#     beamline_element = BeamlineElement(optical_element=optical_element,
#                                        coordinates=ElementCoordinates(p=20.000000, q=0.000000,
#                                                                       angle_radial=numpy.radians(0.000000),
#                                                                       angle_azimuthal=numpy.radians(0.000000)))
#     propagation_elements.add_beamline_element(beamline_element)
#     propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
#     # self.set_additional_parameters(propagation_parameters)
#     #
#     propagation_parameters.set_additional_parameters('magnification_x', 0.5)
#     propagation_parameters.set_additional_parameters('magnification_N', 1.0)
#     #
#     propagator = PropagationManager.Instance()
#     try:
#         propagator.add_propagator(Integral1D())
#     except:
#         pass
#     output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
#                                                  handler_name='INTEGRAL_1D')
#
#     #
#     # ---- plots -----
#     #
#     # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 3')
#
#     ##########  OPTICAL ELEMENT NUMBER 4 ##########
#
#     input_wavefront = output_wavefront.duplicate()
#     from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D
#
#     optical_element = WOIdealLens1D(name='', focal_length=28.200000)
#
#     # no drift in this element
#     output_wavefront = optical_element.applyOpticalElement(input_wavefront)
#
#     #
#     # ---- plots -----
#     #
#     # plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 4')
#
#     ##########  OPTICAL ELEMENT NUMBER 5 ##########
#
#     input_wavefront = output_wavefront.duplicate()
#     from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D
#
#     optical_element = WOScreen1D()
#
#     # drift_after 20 m
#     #
#     # propagating
#     #
#     #
#     propagation_elements = PropagationElements()
#     beamline_element = BeamlineElement(optical_element=optical_element,
#                                        coordinates=ElementCoordinates(p=0.000000, q=q,
#                                                                       angle_radial=numpy.radians(0.000000),
#                                                                       angle_azimuthal=numpy.radians(0.000000)))
#     propagation_elements.add_beamline_element(beamline_element)
#     propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
#     # self.set_additional_parameters(propagation_parameters)
#     #
#     propagation_parameters.set_additional_parameters('magnification_x', 3.0)
#     #
#     propagator = PropagationManager.Instance()
#     try:
#         propagator.add_propagator(FresnelZoom1D())
#     except:
#         pass
#     output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
#                                                  handler_name='FRESNEL_ZOOM_1D')
#
#     return output_wavefront

def run_wofry_h(plot_from=0,p7=24.000000,p8=37.000000):
    DISTANCE = []
    FWHM = []
    I0 = []
    ##########  SOURCE ##########

    #
    # create output_wavefront
    #
    #
    output_wavefront = GenericWavefront1D.initialize_wavefront_from_range(x_min=-0.00012, x_max=0.00012,
                                                                          number_of_points=1000)
    output_wavefront.set_photon_energy(10000)
    output_wavefront.set_gaussian_hermite_mode(sigma_x=3.03783e-05, amplitude=1, mode_x=0, shift=0, beta=0.0922395)

    if plot_from <= 0: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='SOURCE')

    FWHM.append(1e6*get_wavefront_intensity_fwhm(output_wavefront))
    I0.append(get_wavefront_intensity_I0(output_wavefront))
    DISTANCE.append(0.0)
    ##########  OPTICAL SYSTEM ##########

    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 35 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=35.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 8.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    #
    # ---- plots -----
    #
    if plot_from <= 1: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 1')

    FWHM.append(1e6*get_wavefront_intensity_fwhm(output_wavefront))
    I0.append(get_wavefront_intensity_I0(output_wavefront))
    DISTANCE.append(35)

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *
    boundary_shape = Rectangle(-1.25e-05, 1.25e-05, -1.25e-05, 1.25e-05)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    FWHM.append(1e6*get_wavefront_intensity_fwhm(output_wavefront))
    I0.append(get_wavefront_intensity_I0(output_wavefront))
    DISTANCE.append(35.01)
    #
    # ---- plots -----
    #
    if plot_from <= 2: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 2')

    ##########  OPTICAL ELEMENT NUMBER 3 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 20 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=20.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 0.5)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='INTEGRAL_1D')

    FWHM.append(1e6*get_wavefront_intensity_fwhm(output_wavefront))
    I0.append(get_wavefront_intensity_I0(output_wavefront))
    DISTANCE.append(55)
    #
    # ---- plots -----
    #
    if plot_from <= 3: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 3')

    ##########  OPTICAL ELEMENT NUMBER 4 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='', focal_length=28.200000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    #
    # ---- plots -----
    #
    if plot_from <= 4: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 4')

    FWHM.append(1e6*get_wavefront_intensity_fwhm(output_wavefront))
    I0.append(get_wavefront_intensity_I0(output_wavefront))
    DISTANCE.append(55.01)

    ##########  OPTICAL ELEMENT NUMBER 5 ##########

    input_wavefront = output_wavefront.duplicate()
    # from syned.beamline.shape import *
    boundary_shape = Rectangle(-0.0005, 0.0005, -0.0005, 0.0005)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # drift_before 99 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=99.000000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 2.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    FWHM.append(1e6*get_wavefront_intensity_fwhm(output_wavefront))
    I0.append(get_wavefront_intensity_I0(output_wavefront))
    DISTANCE.append(164)
    #
    # ---- plots -----
    #
    if plot_from <= 5: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 5')

    ##########  OPTICAL ELEMENT NUMBER 6 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.lens import WOIdealLens1D

    optical_element = WOIdealLens1D(name='', focal_length=39.700000)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    FWHM.append(1e6*get_wavefront_intensity_fwhm(output_wavefront))
    I0.append(get_wavefront_intensity_I0(output_wavefront))
    DISTANCE.append(164.01)

    #
    # ---- plots -----
    #
    if plot_from <= 6: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 6')

    ##########  OPTICAL ELEMENT NUMBER 7 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 24 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=p7, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    #
    # ---- plots -----
    #
    if plot_from <= 7: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 7')

    ##########  OPTICAL ELEMENT NUMBER 8 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 37 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=p8, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(FresnelZoom1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='FRESNEL_ZOOM_1D')

    #
    # ---- plots -----
    #
    if plot_from <= 8: plot(output_wavefront.get_abscissas(), output_wavefront.get_intensity(), title='OPTICAL ELEMENT NR 8')

    return input_wavefront, FWHM, I0, DISTANCE # oe.7

if __name__ == "__main__":


    # pscan = numpy.linspace(10.0, 300.0, 10)
    # FWHM = numpy.zeros_like(pscan)
    # I0 = numpy.zeros_like(pscan)
    #
    # for i in range(pscan.size):
    #     wf = first_lens_focus(q=pscan[i])
    #
    #     FWHM[i] = 1e6*get_wavefront_intensity_fwhm(wf)
    #     I0[i] = get_wavefront_intensity_I0(wf)
    #     txt = "i: %d p: %g  FWHM: %g um; I0: %g" % (i, pscan[i], FWHM[i], I0[i])
    #     print(txt)
    #     # plot(1e6*wf.get_abscissas(),wf.get_intensity(),title = txt)


    p7 = 24.000000
    p8 = 37.000000
    pscan = numpy.linspace(0.0, 75.0, 10)
    # FWHM = numpy.zeros_like(pscan)
    # I0 = numpy.zeros_like(pscan)
    wf, FWHM, I0, DISTANCE = run_wofry_h(plot_from=18, p7=p7, p8=p8)
    for i in range(pscan.size):
        p7i = p7+pscan[i]
        p8i = (p8-pscan[i])
        wf, _, _, _ = run_wofry_h(plot_from=18,p7=p7i,p8=p8i)

        FWHM.append(1e6*get_wavefront_intensity_fwhm(wf))
        I0.append(get_wavefront_intensity_I0(wf))
        DISTANCE.append(164.0+p7+pscan[i])
        txt = "p7: %g  p8:%g FWHM: %g um; I0: %g" % (p7i, p8i, FWHM[-1], I0[-1])
        print(txt)
        # plot(1e6*wf.get_abscissas(),wf.get_intensity(),title = txt)


    plot(DISTANCE, FWHM)
    plot(DISTANCE, I0)
