from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter

from srxraylib.plot.gol import plot_image, plot

import h5py

from vortx_propagate import AFpropagated, W_at_x2x2, propagate

#
########################################################################################################################
#
#
#
# def W_at_x2x2(af,index_x2=None,index_y2=None,index_max=0):
#
#
#     if index_x2 is None:
#         index_x2 = int(af.x_coordinates().size / 2)
#
#     if index_y2 is None:
#         index_x2 = int(af.y_coordinates().size / 2)
#
#     print("Using indices: ",index_x2,index_y2," out of ",af.x_coordinates().size,af.y_coordinates().size)
#
#     # for i in range(index_max+1):
#     #     print(">>>>>>>>>>>>",i,af.eigenvalue(i),af.mode(i).shape)
#
#     for i in range(index_max+1):
#         # print(">>>>>>> i, index_max",i,index_max)
#         mi = af.mode(i)
#         evi = af.eigenvalue(i)
#
#         if i == 0:
#             # print(">>>>>>>>>>> assigning W")
#             W  = numpy.conj(mi) * evi * mi[index_x2,index_y2]
#         else:
#             W += evi * numpy.conj(mi) * mi[index_x2,index_y2]
#
#     return W
#
# def propagate(af,distance=30.0,index_max=9,zoom=(2.0,6.0)):
#
#
#     from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
#     from syned.beamline.beamline_element import BeamlineElement
#     from syned.beamline.element_coordinates import ElementCoordinates
#     from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
#
#     from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
#     from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen
#
#     afp = AFpropagated()
#
#
#
#     for i in range(index_max+1):
#         mi = af.mode(i)
#         evi = af.eigenvalue(i)
#
#         print("propagating mode index",i,evi,mi.shape)
#
#
#         input_wavefront = GenericWavefront2D.initialize_wavefront_from_arrays(x_array=af.x_coordinates(),
#                                                                               y_array=af.y_coordinates(),
#                                                                               z_array=mi,
#                                                                               )
#
#         input_wavefront.set_photon_energy(17226.0)
#
#
#         optical_element = WOScreen()
#
#
#         propagation_elements = PropagationElements()
#         beamline_element = BeamlineElement(optical_element=optical_element,
#                                         coordinates=ElementCoordinates(p=distance,q=0.000000,
#                                         angle_radial=numpy.radians(0.000000),angle_azimuthal=numpy.radians(0.000000)))
#         propagation_elements.add_beamline_element(beamline_element)
#         propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),propagation_elements=propagation_elements)
#
#         propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
#         propagation_parameters.set_additional_parameters('magnification_x', zoom[0])
#         propagation_parameters.set_additional_parameters('magnification_y', zoom[1])
#
#         propagator = PropagationManager.Instance()
#         try:
#             propagator.add_propagator(FresnelZoomXY2D())
#         except:
#             pass
#
#         output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,handler_name='FRESNEL_ZOOM_XY_2D')
#
#
#         # plot_image(output_wavefront.get_intensity())
#
#         afp.add_mode(output_wavefront.get_complex_amplitude(),evi)
#         if i == 0:
#             afp.x = output_wavefront.get_coordinate_x()
#             afp.y = output_wavefront.get_coordinate_y()
#
#     # for i in range(index_max+1):
#     #     print(">>>>",afp.mode(i).shape,afp.x_coordinates().shape)
#     #     plot_image(afp.get_intensity(i),1e6*afp.x_coordinates(),1e6*afp.y_coordinates(),title="up to i=%d"%i)
#     #
#     # print(">>>>>>>>>>> number of modes in afp: ",afp.number_modes())
#     return afp

if __name__ == "__main__":



    filename_ebs="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"



    point = "C5"

    distance = 30.0
    index_max = 1099
    zoom = (6.0,16.0)

    if point == "A":
        coordinate_x = 0.0
        coordinate_y = 0.0
    elif point == "B":
        coordinate_x = 57.0e-6
        coordinate_y = 148.2e-6
    elif point == "C":
        coordinate_x = 125.47e-6
        coordinate_y = 312.78e-6
    elif point == "D":
        coordinate_x = -7e-6
        coordinate_y = 31e-6
    elif point == "C1":
        distance = 1.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        index_max = 99
        zoom = (1.0,1.0)
    elif point == "C5":
        distance = 5.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        index_max = 999 # 9 # 99 # 19
        zoom = (1.0,1.5)
    elif point == "C100":
        distance = 100.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        index_max = 99
        zoom = (10.0,40.0)
    else:
        raise Exception("Point not found!")



    af  = CompactAFReader.initialize_from_file(filename_ebs)


    SCAN = numpy.linspace(-1,1,9)





    #
    # propagate
    #
    for i in range(SCAN.size):

        print("\n\n######################### propagating point %d of %d\n\n"%(i,SCAN.size-1))
        afp = propagate(af,distance=distance+SCAN[i],index_max=index_max,zoom=zoom)


        #
        # get indices
        #

        print("X: start, step, points",afp.x_coordinates()[0],afp.x_coordinates()[1] - afp.x_coordinates()[0],afp.x_coordinates().size)
        print("Y: start, step, points",afp.y_coordinates()[0],afp.y_coordinates()[1] - afp.y_coordinates()[0],afp.y_coordinates().size)
        step_x = afp.x_coordinates()[1] - afp.x_coordinates()[0]
        step_y = afp.y_coordinates()[1] - afp.y_coordinates()[0]
        origin_x = afp.x_coordinates()[0]
        origin_y = afp.y_coordinates()[0]

        index_x2 = int((coordinate_x - origin_x) / step_x)
        index_y2 = int((coordinate_y - origin_y) / step_y)


        print("Using indices: ",index_x2,index_y2," out of ",afp.x_coordinates().size,afp.y_coordinates().size,
              "ratio: ",index_x2/afp.x_coordinates().size,index_y2/afp.y_coordinates().size)


        if i == 0:
            # output file initialization
            h5file = "vx_id16a_%s_propagated_neighbour_mode%04d.h5"%(point,index_max)
            h5w = H5SimpleWriter.initialize_file(h5file,creator="vortx_propagate_neighbour.py")
            h5w.add_key("r2_indices",[index_x2,index_y2])
            h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]])
            h5w.add_key("SCAN",SCAN)
            h5w.add_key("DISTANCES",numpy.array(SCAN)+distance)


        tmp = W_at_x2x2(afp,index_x2=index_x2,index_y2=index_y2,index_max=index_max)

        #




        x_ebs = numpy.arange(afp.number_modes())
        y_ebs = numpy.abs(afp.occupation_array())


        print("Calculating mode index %d"%index_max)
        tmp = W_at_x2x2(afp,index_x2=index_x2,index_y2=index_y2,index_max=index_max)

        h5w.create_entry("uptomode%04d_%04d"%(index_max,i),nx_default="Wcomplex")
        h5w.add_key("r2_indices",[index_x2,index_y2], entry_name="uptomode%04d_%04d"%(index_max,i))
        h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]], entry_name="uptomode%04d_%04d"%(index_max,i))

        h5w.add_image(tmp,1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d_%04d"%(index_max,i),
                     image_name="Wcomplex",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(numpy.absolute(tmp),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d_%04d"%(index_max,i),
                     image_name="Wamplitude",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(numpy.angle(tmp),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d_%04d"%(index_max,i),
                     image_name="Wphase",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(afp.get_intensity(index_max),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d_%04d"%(index_max,i),
                     image_name="SpectralDensity",title_x="X [mm]",title_y="Y [mm]")

    print("File written to disk: ",h5file)