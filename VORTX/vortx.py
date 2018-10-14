from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter
#
########################################################################################################################
#


def W_at_x2x2(af,index_x2=None,index_y2=None,index_max=0):

    if index_x2 is None:
        index_x2 = int(af.x_coordinates().size / 2)

    if index_y2 is None:
        index_x2 = int(af.y_coordinates().size / 2)

    print("Using indices: ",index_x2,index_y2," out of ",af.x_coordinates().size,af.y_coordinates().size)

    for i in range(index_max+1):
        mi = af.mode(i)
        evi = af.eigenvalue(i)

        if i == 0:
            W  = numpy.conj(mi) * evi * mi[index_x2,index_y2]
        else:
            W += evi * numpy.conj(mi) * mi[index_x2,index_y2]

    return W

if __name__ == "__main__":

    filename_ebs="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"


    af  = CompactAFReader.initialize_from_file(filename_ebs)


    h5file = "vx_id16a_C.h5"
    index_x2=573
    index_y2=200

    #
    # h5file = "vx_id16a_B.h5"
    # index_x2=535
    # index_y2=183
    #

    # h5file = "vx_id16a_A.h5"
    # index_x2=1007//2
    # index_y2=335//2


    if index_x2 is None:
        index_x2 = int(af.x_coordinates().size / 2)

    if index_y2 is None:
        index_y2 = int(af.y_coordinates().size / 2)

    print("Unsing indices: ",index_x2,index_y2," out of ",af.x_coordinates().size,af.y_coordinates().size)


    #


    h5w = H5SimpleWriter.initialize_file(h5file,creator="vortx.py")
    # h5_initialize(h5file,creator="vortx.py")
    h5w.add_key("r2_indices",[index_x2,index_y2])
    h5w.add_key("r2",[af.x_coordinates()[index_x2],af.y_coordinates()[index_y2]])


    x_ebs = numpy.arange(af.number_modes())
    y_ebs = numpy.abs(af.occupation_array())

    # plot(x_ebs,y_ebs)

    # h5w.aadd_dataset()



    # t = numpy.concatenate((numpy.arange(10),numpy.arange(10,100,10),numpy.arange(100,1001,100)))
    t = numpy.array( (0,1,2,3,4,5,6,7,8,9,
                      19,29,39,49,59,69,79,89,99,
                      199,299,399,499,599,699,799,899,999,
                      1099) )

    print(t)

    for index_max in t:
        print("Calculating mode %d"%index_max)
        # tmp = W_at_x2x2(af,index_x2=585,index_y2=170,index_max=index_max)
        tmp = W_at_x2x2(af,index_x2=index_x2,index_y2=index_y2,index_max=index_max)


        h5w.create_entry("uptomode%04d"%index_max,nx_default="SpectralDensity")
        h5w.add_key("r2_indices",[index_x2,index_y2], entry_name="uptomode%04d"%index_max)
        h5w.add_key("r2",[af.x_coordinates()[index_x2],af.y_coordinates()[index_y2]], entry_name="uptomode%04d"%index_max)

        h5w.add_image(tmp,1e3*af.x_coordinates(),1e3*af.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="Wcomplex",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(numpy.absolute(tmp),1e3*af.x_coordinates(),1e3*af.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="Wamplitude",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(af.intensity_from_modes(index_max),1e3*af.x_coordinates(),1e3*af.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="SpectralDensity",title_x="X [mm]",title_y="Y [mm]")

        h5w.add_image(numpy.angle(tmp),1e3*af.x_coordinates(),1e3*af.y_coordinates(),
                     entry_name="uptomode%04d"%index_max,
                     image_name="Wphase",title_x="X [mm]",title_y="Y [mm]")
