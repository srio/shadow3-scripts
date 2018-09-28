from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter

from srxraylib.plot.gol import plot_image, plot

import h5py

from vortx_propagate import AFpropagated, W_at_x2x2, propagate



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
        index_max = 9 # 99 # 19
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