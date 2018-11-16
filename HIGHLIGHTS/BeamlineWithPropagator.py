
#
# Extend the a syned beamline (filled with wofry elements) with propagator parameters and makes propagation of
# wavefronts and af (comsyl)
#


import numpy

from wofry.propagator.propagator import PropagationElements, PropagationParameters
from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D


from syned.beamline.beamline import Beamline
from syned.storage_ring.light_source import LightSource


class BeamlineWithPropagator(Beamline):
    def __init__(self,light_source=LightSource(), beamline_elements_list=None, propagator_specific_parameteres=None):
        # self._propagator_handlers = []
        if propagator_specific_parameteres is None:
            self._propagator_specific_parameteres = []
        else:
            self._propagator_specific_parameteres = propagator_specific_parameteres
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)

    def get_propagator_specific_parameteres(self):
        return self._propagator_specific_parameteres

    def set_propagator_specific_parameteres(self,propagator_specific_parameteres):
        self._propagator_specific_parameteres = propagator_specific_parameteres

    def append_propagator_specific_parameteres(self, value):
            self._propagator_specific_parameteres.append(value)

    def get_propagator_specific_parameteres_number(self):
        return len(self._propagator_specific_parameteres)

    def get_propagator_specific_parameteres_at(self, index):
        if index >= len(self._propagator_specific_parameteres):
            raise IndexError("Index " + str(index) + " out of bounds")

        return self._propagator_specific_parameteres[index]

    def append_beamline_element_and_propagator_specific_parameteres(self, element, specific):
        self.append_beamline_element(element)
        self.append_propagator_specific_parameteres(specific)


    # overwrite these methods
    def duplicate(self):
        beamline_elements_list = []
        propagator_specific_parameteres_list = []
        self._check_consistency()
        for beamline_element,specific_parameter in zip(self._beamline_elements_list,
                            self._propagator_specific_parameteres):
            beamline_elements_list.append(beamline_element)
            propagator_specific_parameteres_list.append(specific_parameter)

        return BeamlineWithPropagator(light_source=self._light_source,
                        beamline_elements_list=beamline_elements_list,
                        propagator_specific_parameteres=propagator_specific_parameteres_list)

    def _check_consistency(self):
        if len(self._beamline_elements_list) != len(self._propagator_specific_parameteres):
            raise Exception("Bad number of specific parameters")


    def propagate(self,wofry_wavefront,mypropagator,return_wavefront_for_every_element=False):

        if return_wavefront_for_every_element:
            w_out = wofry_wavefront.duplicate()
            output_wavefronts = []
            output_wavefronts.append(wofry_wavefront.duplicate())
            for i in range(self.get_beamline_elements_number()):
                w_in = w_out.duplicate()

                #
                # propagating
                #
                #
                propagation_elements = PropagationElements()
                propagation_elements.add_beamline_element(self.get_beamline_element_at(i),self.get_propagator_specific_parameteres_at(i))

                propagation_parameters = PropagationParameters(wavefront=w_in,
                                                               propagation_elements = propagation_elements)

                w_out = mypropagator.do_propagation(propagation_parameters=propagation_parameters,
                                                        handler_name='FRESNEL_ZOOM_XY_2D')

                output_wavefronts.append(w_out)

            return output_wavefronts
        else:
            propagation_elements = PropagationElements()
            for i in range(self.get_beamline_elements_number()):
                propagation_elements.add_beamline_element(self.get_beamline_element_at(i),
                                                          self.get_propagator_specific_parameteres_at(i))

            propagation_parameters = PropagationParameters(wavefront=wofry_wavefront.duplicate(),
                                                               propagation_elements = propagation_elements)

            w_out = mypropagator.do_propagation(propagation_parameters=propagation_parameters,
                                                    handler_name='FRESNEL_ZOOM_XY_2D')

            return w_out


    @classmethod
    def _numpy_wavefront_file_to_wofry_wavefront(cls,file_in):
        return cls._numpy_wavefront_to_wofry_wavefront(numpy.load(file_in))

    @classmethod
    def _numpy_wavefront_to_wofry_wavefront(cls,dict_in):
        e_field = dict_in["e_field"]
        coordinates = dict_in["coordinates"]
        energies = dict_in["energies"]

        x = numpy.linspace(coordinates[0],coordinates[1],e_field.shape[1])
        y = numpy.linspace(coordinates[2],coordinates[3],e_field.shape[2])
        wofry_wf_in = GenericWavefront2D.initialize_wavefront_from_arrays(x,y,e_field[0,:,:,0].copy())
        wofry_wf_in.set_photon_energy(energies[0])
        return wofry_wf_in

    @classmethod
    def _wofry_wavefront_to_numpy_wavefront(cls,wofry_wf_out):

        ca = wofry_wf_out.get_complex_amplitude()

        e_field = numpy.zeros((1,ca.shape[0],ca.shape[1],2),dtype=complex)
        e_field[0,:,:,0] = ca

        coordinates = numpy.zeros(4)
        coordinates[0] = wofry_wf_out.get_coordinate_x()[0]
        coordinates[1] = wofry_wf_out.get_coordinate_x()[-1]
        coordinates[2] = wofry_wf_out.get_coordinate_y()[0]
        coordinates[3] = wofry_wf_out.get_coordinate_y()[-1]

        energies = numpy.array([wofry_wf_out.get_photon_energy()])

        return {"e_field":e_field,"coordinates":coordinates,"energies":energies}


    def propagate_numpy_wavefront(self,dict_in,mypropagator):

        wofry_wf_in = self._numpy_wavefront_to_wofry_wavefront(dict_in)

        wofry_wf_out = self.propagate(wofry_wf_in,mypropagator,return_wavefront_for_every_element=False)


        return self._wofry_wavefront_to_numpy_wavefront(wofry_wf_out)


    def propagate_numpy_wavefront_from_files(self,filename_in,filename_out,mypropagator):

        dict_content = numpy.load(filename_in)

        dict_content_out = self.propagate_numpy_wavefront(dict_content,mypropagator)

        numpy.savez(filename_out,
                 e_field     = dict_content_out["e_field"],
                 coordinates = dict_content_out["coordinates"],
                 energies    = dict_content_out["energies"])


        return dict_content_out

    def propagate_af(self,af):
        raise NotImplementedError

def create_beamline_wofry(load_from_file=None, slit_width=5e-6,slit_height=5e-6,source_offset=0.0):

    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen
    from wofry.beamline.optical_elements.ideal_elements.lens import WOIdealLens
    from wofry.beamline.optical_elements.absorbers.slit import WOSlit
    from syned.beamline.shape import Rectangle
    from wofry.propagator.propagator import PropagationElements, PropagationParameters
    # from syned.beamline.beamline import Beamline

    from syned.storage_ring.empty_light_source import EmptyLightSource

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


    #
    # beamline
    #

    # beamline1 = Beamline(beamline_elements_list=[])
    beamline1 = BeamlineWithPropagator()
    beamline1.set_light_source(EmptyLightSource())


    for k in range(len(elements)):
        beamline1.append_beamline_element_and_propagator_specific_parameteres(elements[k],specific[k])

    return beamline1


def test_propagate_wavefront():
    from srxraylib.plot.gol import plot_image, plot

    from wofry.propagator.propagator import PropagationManager
    from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
    from comsyl.utils.Logger import log

    mypropagator = PropagationManager.Instance()
    try:
        mypropagator.add_propagator(FresnelZoomXY2D())
    except:
        log("May be you already initialized propagator and stored FresnelZoomXY2D")


    # BEAMLINE = pickle.load(open("/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/BEAMLINE.p","rb"))
    BEAMLINE = create_beamline_wofry(slit_width=25e-6,slit_height=25e-6,source_offset=1.0)

    print(BEAMLINE.info())

    print(BEAMLINE.get_propagator_specific_parameteres())

    BEAMLINE._check_consistency()
    BEAMLINE.to_json("bl.json")

    file_in =  "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/tmp-working/tmp0_hib3-3302_IN.npz"
    file_out = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/MARK/tmp-working/tmp0_hib3-3302_OUT2.npz"
    wf_out = BEAMLINE.propagate_numpy_wavefront_from_files(file_in,file_out,mypropagator)


    wf_in = BeamlineWithPropagator._numpy_wavefront_file_to_wofry_wavefront(file_in)
    plot_image(wf_in.get_intensity(),
               wf_in.get_coordinate_x()*1e6,
               wf_in.get_coordinate_y()*1e6,title="IN",)

    wf_out = BeamlineWithPropagator._numpy_wavefront_file_to_wofry_wavefront(file_out)
    plot_image(wf_out.get_intensity(),
               wf_out.get_coordinate_x()*1e6,
               wf_out.get_coordinate_y()*1e6,title="OUT")


def test_propagate_af():
    from srxraylib.plot.gol import plot_image, plot

    from wofry.propagator.propagator import PropagationManager
    from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
    from comsyl.utils.Logger import log
    from comsyl.autocorrelation.AutocorrelationFunction import AutocorrelationFunction

    mypropagator = PropagationManager.Instance()
    try:
        mypropagator.add_propagator(FresnelZoomXY2D())
    except:
        log("May be you already initialized propagator and stored FresnelZoomXY2D")


    # BEAMLINE = pickle.load(open("/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/BEAMLINE.p","rb"))
    BEAMLINE = create_beamline_wofry(slit_width=25e-6,slit_height=25e-6,source_offset=1.0)

    print(BEAMLINE.info())

    print(BEAMLINE.get_propagator_specific_parameteres())

    BEAMLINE._check_consistency()
    BEAMLINE.to_json("bl.json")

    filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    af_name = filename.split("/")[-1].replace(".npz", "")
    autocorrelation_function = AutocorrelationFunction.load(filename)


    af_out = BEAMLINE.propagate_af(autocorrelation_function)

    # af_out.save("")


if __name__ == "__main__":


    # test_propagate_wavefront()

    test_propagate_af()
