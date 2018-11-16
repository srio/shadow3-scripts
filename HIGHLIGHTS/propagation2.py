#
# runs using "WofrySuperBeamline" with wofry wavefronts
#


from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
# from CompactAFReader import CompactAFReader

import numpy


from srxraylib.plot.gol import plot_image, plot

from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
from wofry.propagator.propagator import PropagationManager


from wofry.propagator.propagator import PropagationElements, PropagationParameters
from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen
from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens
from wofry.beamline.optical_elements.absorbers.slit import WOSlit

from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D

from syned.beamline.shape import Rectangle

class WofrySuperBeamline(object):
    def __init__(self):
        self._beamline_elements = []
        self._propagator_handlers = []
        self._propagator_specific_parameteres = []

        self._initialize_propagator()

    def add_element(self,
                    beamline_element=BeamlineElement(),
                    propagator_handler="FRESNEL_ZOOM_XY_2D",
                    propagator_specific_parameters={'shift_half_pixel':1,'magnification_x':1.0,'magnification_y':1.0},
                    ):
        self._beamline_elements.append(beamline_element)
        self._propagator_handlers.append(propagator_handler)
        self._propagator_specific_parameteres.append(propagator_specific_parameters)


    def number_of_elements(self):
        return len(self._beamline_elements)

    def _initialize_propagator(self):

        self._propagator = PropagationManager.Instance()
        try:
            self._propagator.add_propagator(FresnelZoomXY2D())
        except:
            print("May be you alreay initialized propagator and stored FresnelZoomXY2D")

    def _info(self):

        for i in range(self.number_of_elements()):
            print(">>>",self._beamline_elements[i])

    def propagate(self,wofry_wavefront):

        w_out = wofry_wavefront.duplicate()
        output_wavefronts = []

        for i in range(self.number_of_elements()):
            w_in = w_out.duplicate()

            #
            # propagating
            #
            #
            propagation_elements = PropagationElements()
            propagation_elements.add_beamline_element(self._beamline_elements[i])

            propagation_parameters = PropagationParameters(wavefront=w_in,propagation_elements = propagation_elements)

            propagation_parameters.set_additional_parameters('shift_half_pixel',(self._propagator_specific_parameteres[i])["shift_half_pixel"])
            propagation_parameters.set_additional_parameters('magnification_x', (self._propagator_specific_parameteres[i])["magnification_x"])
            propagation_parameters.set_additional_parameters('magnification_y', (self._propagator_specific_parameteres[i])["magnification_y"])

            w_out = self._propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                    handler_name=self._propagator_handlers[i])

            output_wavefronts.append(w_out)



        return output_wavefronts


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

def create_wofry_elements(slit_width=5e-6,slit_height=5e-6):

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

    return [beamline_element1,beamline_element2,beamline_element3,beamline_element4],\
           ['FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D','FRESNEL_ZOOM_XY_2D'],\
           [{'shift_half_pixel':1,'magnification_x':8.0,'magnification_y':8.0},
            {'shift_half_pixel':1,'magnification_x':1.0,'magnification_y':1.0},
            {'shift_half_pixel':1,'magnification_x':0.5,'magnification_y':0.35},
            {'shift_half_pixel':1,'magnification_x':0.2,'magnification_y':0.2}]


if __name__ == "__main__":


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
            e,h,s = create_wofry_elements(slit_width=slit_width,slit_height=slit_height)

            # plot_image(wofry_wavefront.get_intensity(),1e6*wofry_wavefront.get_coordinate_x(),1e6*wofry_wavefront.get_coordinate_y(),
            #            title="source")


            BEAMLINE = None
            BEAMLINE = WofrySuperBeamline()

            for k in range(len(e)):
                BEAMLINE.add_element(beamline_element=e[k],propagator_handler=h[k],propagator_specific_parameters=s[k])

            # BEAMLINE._info()

            wf_list = BEAMLINE.propagate(wofry_wavefront)



            intensity_accumulated += wf_list[-1].get_integrated_intensity()

            # for ll in range(4):
            #     plot_image(wf_list[ll].get_intensity(),1e6*wf_list[ll].get_coordinate_x(),1e6*wf_list[ll].get_coordinate_y(),
            #                title="slit H:%d x V:%d ; mode %d propagated to final oe"%(slitH[j],slitV[j],i)+">>ll% d"%ll)

            ratio = intensity_accumulated/intensity_full
            print("slit H:%d, V:%d, up to mode: %d intensity(full): %g intensity(cut): %g ratio: %g "%(
                slitH[j],slitV[j],i,intensity_full,intensity_accumulated,ratio,))


            print("\n\n\n\n\n\n")
        f = open("propagation.txt",'a')
        f.write('"%d x %d" %f \n'%(slitV[j],slitH[j],ratio))
        f.close
