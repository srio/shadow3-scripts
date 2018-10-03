from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter

from srxraylib.plot.gol import plot_image, plot

import h5py

#
########################################################################################################################
#

class AFpropagated(object):
    def __init__(self):
        self.x = None
        self.y = None
        self.modes = []
        self.eigenvalues = []

    def add_mode(self,mode,eigenvalue):
        self.modes.append(mode)
        self.eigenvalues.append(eigenvalue)

    def x_coordinates(self):
        return self.x.copy()

    def y_coordinates(self):
        return self.y.copy()

    def mode(self,index):
        return self.modes[index].copy()

    def eigenvalue(self,index):
        return self.eigenvalues[index]

    def number_modes(self):
        return len(self.eigenvalues)

    def occupation_array(self):
        return numpy.array(self.eigenvalues).real

    def save_h5file(self,h5file):

        # h5file = "vx_id16a_A_propagated.h5"

        h5w = H5SimpleWriter.initialize_file(h5file,creator="vortx_propagate.py")
        # h5w.add_key("r2_indices",[index_x2,index_y2])
        # h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]])


        x_ebs = numpy.arange(self.number_modes())
        y_ebs = numpy.abs(self.occupation_array())

        # plot(x_ebs,y_ebs)


        for index in range(self.number_modes()):

            h5w.create_entry("mode%04d"%index,nx_default="Wintensity")
            # h5w.add_key("r2_indices",[index_x2,index_y2], entry_name="uptomode%04d"%index_max)
            # h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]], entry_name="uptomode%04d"%index_max)

            h5w.add_image(self.mode(index),1e3*self.x_coordinates(),1e3*self.y_coordinates(),
                         entry_name="mode%04d"%index,
                         image_name="Wcomplex",title_x="X [mm]",title_y="Y [mm]")

            h5w.add_image(self.get_intensity(index),1e3*self.x_coordinates(),1e3*self.y_coordinates(),
                         entry_name="mode%04d"%index,
                         image_name="Wintensity",title_x="X [mm]",title_y="Y [mm]")


        h5w.create_entry("eigenvalues",nx_default="occupation")
        h5w.add_dataset(x_ebs,y_ebs,entry_name="eigenvalues",dataset_name="occupation")
        # h5w.add_image(self.get_intensity(5),1e3*self.x_coordinates(),1e3*self.y_coordinates(),
        #              entry_name="eigenvalues",
        #              image_name="occupation",title_x="X [mm]",title_y="Y [mm]")



    @classmethod
    def load_h5file(cls,h5file):

        f = h5py.File(h5file,'r')

        eigenvalues = f["eigenvalues/occupation/y"].value

        oo = AFpropagated()

        for i in range(eigenvalues.size):
            oo.add_mode(f["mode%04d/Wcomplex/image_data"%i].value.T,eigenvalues[i])
            if i == 0:
                oo.x = f["mode0000/Wcomplex/axis_x"].value
                oo.y = f["mode0000/Wcomplex/axis_y"].value

        f.close()

        return oo

    def info(self):
        print(">> number of modes: ",self.number_modes(),len(self.modes))
        print(">> shape of mode: ",self.mode(0).shape)
        print(">> shape of x,y: ",self.x.shape,self.y.shape)
        print(">> max,min of 0 mode: ",self.mode(0).min(),self.mode(0).max())

    # def get_complex_amplitude(self,index_max=None):
    #     if index_max == None:
    #         index_max = self.number_modes() - 1
    #
    #     for i in range(index_max+1):
    #         if i == 0:
    #             tmp = numpy.sqrt(self.eigenvalue(i).real) * self.mode(i)
    #         else:
    #             tmp += numpy.sqrt(self.eigenvalue(i).real) * self.mode(i)
    #     return tmp
    #
    # def get_amplitude(self,index_max=None):
    #     return numpy.abs(self.get_complex_amplitude(index_max=index_max))

    def get_intensity(self,index_max=None):
        if index_max == None:
            index_max = self.number_modes() - 1

        for i in range(index_max+1):
            if i == 0:
                tmp1 = numpy.sqrt(self.eigenvalue(i).real) * self.mode(i)
                tmp = numpy.abs(tmp1)**2
            else:
                tmp1 = numpy.sqrt(self.eigenvalue(i).real) * self.mode(i)
                tmp += numpy.abs(tmp1)**2
        return tmp


def W_at_x2x2(af,index_x2=None,index_y2=None,index_max=0):


    if index_x2 is None:
        index_x2 = int(af.x_coordinates().size / 2)

    if index_y2 is None:
        index_x2 = int(af.y_coordinates().size / 2)

    print("Using indices: ",index_x2,index_y2," out of ",af.x_coordinates().size,af.y_coordinates().size)

    # for i in range(index_max+1):
    #     print(">>>>>>>>>>>>",i,af.eigenvalue(i),af.mode(i).shape)

    for i in range(index_max+1):
        # print(">>>>>>> i, index_max",i,index_max)
        mi = af.mode(i)
        evi = af.eigenvalue(i)

        if i == 0:
            # print(">>>>>>>>>>> assigning W")
            W  = numpy.conj(mi) * evi * mi[index_x2,index_y2]
        else:
            W += evi * numpy.conj(mi) * mi[index_x2,index_y2]

    return W

def propagate(af,distance=30.0,index_max=9,zoom=(2.0,6.0)):


    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    from syned.beamline.beamline_element import BeamlineElement
    from syned.beamline.element_coordinates import ElementCoordinates
    from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D

    from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
    from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen

    afp = AFpropagated()



    for i in range(index_max+1):
        mi = af.mode(i)
        evi = af.eigenvalue(i)

        print("propagating mode index",i,evi,mi.shape)


        input_wavefront = GenericWavefront2D.initialize_wavefront_from_arrays(x_array=af.x_coordinates(),
                                                                              y_array=af.y_coordinates(),
                                                                              z_array=mi,
                                                                              )

        input_wavefront.set_photon_energy(17226.0)


        optical_element = WOScreen()


        propagation_elements = PropagationElements()
        beamline_element = BeamlineElement(optical_element=optical_element,
                                        coordinates=ElementCoordinates(p=distance,q=0.000000,
                                        angle_radial=numpy.radians(0.000000),angle_azimuthal=numpy.radians(0.000000)))
        propagation_elements.add_beamline_element(beamline_element)
        propagation_parameters = PropagationParameters(wavefront=input_wavefront.duplicate(),propagation_elements=propagation_elements)

        propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
        propagation_parameters.set_additional_parameters('magnification_x', zoom[0])
        propagation_parameters.set_additional_parameters('magnification_y', zoom[1])

        propagator = PropagationManager.Instance()
        try:
            propagator.add_propagator(FresnelZoomXY2D())
        except:
            pass

        output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,handler_name='FRESNEL_ZOOM_XY_2D')


        # plot_image(output_wavefront.get_intensity())

        afp.add_mode(output_wavefront.get_complex_amplitude(),evi)
        if i == 0:
            afp.x = output_wavefront.get_coordinate_x()
            afp.y = output_wavefront.get_coordinate_y()

    # for i in range(index_max+1):
    #     print(">>>>",afp.mode(i).shape,afp.x_coordinates().shape)
    #     plot_image(afp.get_intensity(i),1e6*afp.x_coordinates(),1e6*afp.y_coordinates(),title="up to i=%d"%i)
    #
    # print(">>>>>>>>>>> number of modes in afp: ",afp.number_modes())
    return afp


def apply_two_apertures(af,index_max=None,patch_shape='Rectangle',
                center1=[0.0,0.0],width1=[15e-6,15e-6],center2=[40e-6,25e-6],width2=[15e-6,15e-6]):


    from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
    afp = AFpropagated()

    if index_max is None:
        index_max = af.number_modes()

    for i in range(index_max+1):
        mi = af.mode(i)
        evi = af.eigenvalue(i)

        print("applying aperture on mode index",i,evi,mi.shape)


        input_wavefront = GenericWavefront2D.initialize_wavefront_from_arrays(x_array=af.x_coordinates(),
                                                                              y_array=af.y_coordinates(),
                                                                              z_array=mi,
                                                                              )
        input_wavefront.set_photon_energy(17226.0)

        output_wavefront = input_wavefront.duplicate()

        if patch_shape == "Rectangle":
                w1 = output_wavefront.clip_square(center1[0]-0.5*width1[0],center1[0]+0.5*width1[0],
                                            center1[1]-0.5*width1[1],center1[1]+0.5*width1[1],
                                            apply_to_wavefront=False)


                w2 = output_wavefront.clip_square(center2[0]-0.5*width2[0],center2[0]+0.5*width2[0],
                                            center2[1]-0.5*width2[1],center2[1]+0.5*width2[1],
                                            apply_to_wavefront=False)

        elif (patch_shape == "Circle" or patch_shape == "Ellipse"):
                w1 = output_wavefront.clip_circle(width1[0],center1[0],center1[1],
                                            apply_to_wavefront=False)


                w2 = output_wavefront.clip_circle(width2[0],center2[0],center2[1],
                                            apply_to_wavefront=False)

                # plot_image(w1+w2,output_wavefront.get_coordinate_x()*1e6,
                #            output_wavefront.get_coordinate_y()*1e6)
        else:
            raise Exception(NotImplementedError)

        output_wavefront.clip_window(w1+w2)


        afp.add_mode(output_wavefront.get_complex_amplitude(),evi)

        if i == 0:
            afp.x = output_wavefront.get_coordinate_x()
            afp.y = output_wavefront.get_coordinate_y()

    return afp




if __name__ == "__main__":



    filename_ebs="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"



    point = "C100"

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
        zoom = (1.0,1.0)
    elif point == "C5":
        distance = 5.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        zoom = (2.0,2.0)
    elif point == "C100":
        distance = 100.0
        coordinate_x = 125.47e-6 / 30.0 * distance
        coordinate_y = 312.78e-6 / 30.0 * distance
        zoom = (10.0,40.0)
    else:
        raise Exception("Point not found!")



    af  = CompactAFReader.initialize_from_file(filename_ebs)



    # afp = propagate(af,distance=10.0,index_max=index_max,zoom=(2.0,5.0))
    #
    #
    # # afp.save_h5file("id16a_modes_propagated_30m.h5")
    #
    #
    # # print("Loading file: id16a_modes_propagated_30m.h5 ")
    # # afp = AFpropagated.load_h5file("id16a_modes_propagated_30m.h5")
    # # print("Done loading file: id16a_modes_propagated_30m.h5 ")
    # # afp.info()
    # # index_max = afp.number_modes() - 1


    #
    # get indices
    #

    # first propagate a few modes only to check there are no errors
    afp = propagate(af,distance=distance,index_max=5,zoom=zoom)

    h5file = "vx_id16a_%s_propagated.h5"%point

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



    #
    # propagate
    #



    # now propagate all modes
    afp = propagate(af,distance=distance,index_max=index_max,zoom=zoom)


    tmp = W_at_x2x2(afp,index_x2=index_x2,index_y2=index_y2,index_max=index_max)


    #

    h5w = H5SimpleWriter.initialize_file(h5file,creator="vortx_propagate.py")
    # h5_initialize(h5file,creator="vortx.py")
    h5w.add_key("r2_indices",[index_x2,index_y2])
    h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]])


    x_ebs = numpy.arange(afp.number_modes())
    y_ebs = numpy.abs(afp.occupation_array())

    # plot(x_ebs,y_ebs)

    if index_max > 100:
        t = numpy.array( (0,1,2,3,4,5,6,7,8,9,
                          19,29,39,49,59,69,79,89,99,
                          199,299,399,499,599,699,799,899,999,
                          1099) )
    else:
        t = numpy.array( (0,1,2,3,4,5,6,7,8,9,
                          19,29,39,49,59,69,79,89,99, ))


    # t = numpy.array( (0,1,2,3,4,5,6,7,8,9) )

    for index_max in t:
        print("Calculating mode index %d"%index_max)
        tmp = W_at_x2x2(afp,index_x2=index_x2,index_y2=index_y2,index_max=index_max)

        h5w.create_entry("uptomode%04d"%index_max,nx_default="SpectralDensity")
        h5w.add_key("r2_indices",[index_x2,index_y2], entry_name="uptomode%04d"%index_max)
        h5w.add_key("r2",[afp.x_coordinates()[index_x2],afp.y_coordinates()[index_y2]], entry_name="uptomode%04d"%index_max)

        h5w.add_image(tmp,1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="Wcomplex",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(numpy.absolute(tmp),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="Wamplitude",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(numpy.angle(tmp),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="Wphase",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(afp.get_intensity(index_max),1e3*afp.x_coordinates(),1e3*afp.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="SpectralDensity",title_x="X [mm]",title_y="Y [mm]")

    print("File written to disk: ",h5file)