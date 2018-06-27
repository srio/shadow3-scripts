
from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
from syned.util.json_tools import load_from_json_file
from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
from srxraylib.plot.gol import plot_image

from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D



def get_wavefront_of_mode(af,MODE_INDEX):

    #
    # retrieve mode and store it in a generic wavefront
    #
    wf = GenericWavefront2D.initialize_wavefront_from_arrays(
            af.x_coordinates(),af.y_coordinates(), af.mode(MODE_INDEX)  )
    wf.set_photon_energy(af.photon_energy())
    ampl = wf.get_complex_amplitude()
    eigen = af.eigenvalues()
    wf.set_complex_amplitude(ampl * eigen[MODE_INDEX])

    return wf

def propagate_wofry_wavefront_along_beamline(wfr,bl):
    # define propagator to be used
    propagator = PropagationManager.Instance()

    try:
        propagator.add_propagator(FresnelZoomXY2D())
    except:
        pass

    for i,beamline_element in enumerate(bl.get_beamline_elements()): #bl.get_beamline_elements_number()):
        #
        # propagating single element
        #
        print("  >> Propagating element %d of %d"%(i+1,bl.get_beamline_elements_number()))
        propagation_elements = PropagationElements()
        propagation_elements.add_beamline_element(beamline_element)
        propagation_parameters = PropagationParameters(wavefront=wfr.duplicate(),propagation_elements = propagation_elements)

        propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
        propagation_parameters.set_additional_parameters('magnification_x', magnification_x[i])
        propagation_parameters.set_additional_parameters('magnification_y', magnification_y[i])

        wfr = propagator.do_propagation(propagation_parameters=propagation_parameters,handler_name='FRESNEL_ZOOM_XY_2D')

    return wfr

if __name__ == "__main__":

    from srxraylib.plot.gol import plot_image
    

    af = CompactAFReader.initialize_from_file("/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npz")

    # plot_image(wfr.get_intensity(),1e6*wfr.get_coordinate_x(),1e6*wfr.get_coordinate_y(),title="Source",xtitle="X [um]",ytitle="Y [um]")


    bl = load_from_json_file("beamline_id16a.json")


    #
    # propagate mode using Fresnel Zoom propagator
    #
    magnification_x = [5.0, 1,1,0.01,440,1,  1,1,0.00007]
    magnification_y = [10.0,1,1,1,   5,  1,0.5,1,0.00009]

    mode_indices = [0,1,10,100]


    for mode_index in mode_indices:

        print(">>>> propagating mode index: ",mode_index)

        wfr = get_wavefront_of_mode(af,mode_index)

        wfr_final = propagate_wofry_wavefront_along_beamline(wfr,bl)


        plot_image(wfr_final.get_intensity(),1e6*wfr_final.get_coordinate_x(),1e6*wfr_final.get_coordinate_y(),
                   title="Image of mode index: %d"%mode_index,xtitle="X [um]",ytitle="Y [um]")
