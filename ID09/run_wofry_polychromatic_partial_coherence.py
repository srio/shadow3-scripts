


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


#
# SOURCE========================
#


# def run_source(my_mode_index=0):
def run_source(my_mode_index=0,energy=20016.1):

    global coherent_mode_decomposition
    try:
        if my_mode_index == 0: raise Exception()
        tmp = coherent_mode_decomposition
    except:

        ##########  SOURCE ##########

        #
        # create output_wavefront
        #
        #
        from wofryimpl.propagator.util.undulator_coherent_mode_decomposition_1d import \
            UndulatorCoherentModeDecomposition1D
        coherent_mode_decomposition = UndulatorCoherentModeDecomposition1D(
            electron_energy=6,
            electron_current=0.2,
            undulator_period=0.017,
            undulator_nperiods=117.647,
            K=0.09683,
            photon_energy= energy,
            abscissas_interval=0.0001,
            number_of_points=2500,
            distance_to_screen=100,
            scan_direction='V',
            sigmaxx=3.63641e-06,
            sigmaxpxp=1.37498e-06,
            useGSMapproximation=False, )
        # make calculation
        coherent_mode_decomposition_results = coherent_mode_decomposition.calculate()

        mode_index = 0
        output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(mode_index)
    output_wavefront = coherent_mode_decomposition.get_eigenvector_wavefront(my_mode_index)
    return output_wavefront


#
# BEAMLINE========================
#


def run_beamline(output_wavefront):
    ##########  OPTICAL SYSTEM ##########

    ##########  OPTICAL ELEMENT NUMBER 1 ##########

    input_wavefront = output_wavefront.duplicate()
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

    optical_element = WOScreen1D()

    # drift_before 27.066 m
    #
    # propagating
    #
    #
    propagation_elements = PropagationElements()
    beamline_element = BeamlineElement(optical_element=optical_element,
                                       coordinates=ElementCoordinates(p=27.066000, q=0.000000,
                                                                      angle_radial=numpy.radians(0.000000),
                                                                      angle_azimuthal=numpy.radians(0.000000)))
    propagation_elements.add_beamline_element(beamline_element)
    propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
    # self.set_additional_parameters(propagation_parameters)
    #
    propagation_parameters.set_additional_parameters('magnification_x', 20.0)
    propagation_parameters.set_additional_parameters('magnification_N', 1.0)
    #
    propagator = PropagationManager.Instance()
    try:
        propagator.add_propagator(Integral1D())
    except:
        pass
    output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                 handler_name='INTEGRAL_1D')

    ##########  OPTICAL ELEMENT NUMBER 2 ##########

    input_wavefront = output_wavefront.duplicate()
    from syned.beamline.shape import Rectangle
    boundary_shape = Rectangle(-0.0005, 0.0005, -0.0005, 0.0005)
    from wofryimpl.beamline.optical_elements.absorbers.slit import WOSlit1D
    optical_element = WOSlit1D(boundary_shape=boundary_shape)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)

    ##########  OPTICAL ELEMENT NUMBER 3 ##########

    input_wavefront = output_wavefront.duplicate()

    from orangecontrib.esrf.wofry.util.mirror import WOMirror1D

    optical_element = WOMirror1D.create_from_keywords(
        name='',
        shape=0,
        p_focus=44.54,
        q_focus=45.4695,
        grazing_angle_in=0.0025,
        p_distance=17.474,
        q_distance=11.3,
        zoom_factor=2,
        error_flag=1,
        error_file='/home/srio/Oasys/dabam_profile_140461924578000.dat',
        error_file_oversampling_factor=30,
        mirror_length=0,
        mirror_points=0,
        write_profile=0)

    # no drift in this element
    output_wavefront = optical_element.applyOpticalElement(input_wavefront)
    return output_wavefront


#
# MAIN FUNCTION========================
#


# def main():
def main(energy=20016.064):
    from srxraylib.plot.gol import plot, plot_image
    from orangecontrib.esrf.wofry.util.tally import TallyCoherentModes

    tally = TallyCoherentModes()
    for my_mode_index in range(10):
        output_wavefront = run_source(my_mode_index=my_mode_index,energy=energy)
        output_wavefront = run_beamline(output_wavefront)
        tally.append(output_wavefront)

    # tally.plot_cross_spectral_density(show=1, filename="")
    # tally.plot_spectral_density(show=1, filename="")
    # tally.plot_occupation(show=1, filename="")

    tally.save_spectral_density(filename="id09_3mrad_spectral_density.dat")
    tally.save_occupation(filename="id09_3mrad_occupation.dat")


#
# MAIN========================
#


main()

#
# MAIN========================
#

import os
# Energy = numpy.linspace(18000,22000,50)
Energy = numpy.linspace(18500,20500,100)
for energy in Energy:
    main(energy)
    command = "mv id09_3mrad_spectral_density.dat results/id09_3mrad_spectral_density_%4d.dat" % energy
    print(command)
    os.system(command)
    command = "mv id09_3mrad_occupation.dat results/occupation_%4d.dat" % energy
    print(command)
    os.system(command)