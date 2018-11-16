from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
# from CompactAFReader import CompactAFReader

import numpy


from srxraylib.plot.gol import plot_image, plot


# from plot_color import plot_with_transparency_one

import pylab as plt
from matplotlib.colors import Normalize, ListedColormap
import matplotlib.patches as patches


def convert_to_h5(file_from,file_to):
    af = CompactAFReader.initialize_from_file(file_from)
    af.write_h5(file_to)
    print("File written to disk: ",file_to)


if __name__ == "__main__":


    # filename_ebs = "/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"
    # filename_ebs = "/scisoft/data/srio/COMSYL/CALCULATIONS/cs_new_u18_2m_1h_s2.5.h5" # NOT GOOD


    # convert_to_h5("/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz",
    #               "cs_new_u18_2m_1h_s2.5.h5")
    # convert_to_h5("/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy",
    #               "cl_low_beta_u18_2m_1h_s6.5.h5")



    # filename_ebs = "cs_new_u18_2m_1h_s2.5.h5"
    # filename_ebs = "cl_low_beta_u18_2m_1h_s6.5.h5"

    # filename_ebs = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/new_u18_2m_1h_ts_s2.0.npz"
    filename_ebs = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    filename_lb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy" # OK LB
    filename_hb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s2.0.npy"

    #
    # load CSD
    #


    af_ebs  = CompactAFReader.initialize_from_file(filename_ebs)
    cumulated_occupation_ebs = af_ebs.cumulated_occupation_array()
    occupation_ebs = af_ebs.occupation_array()

    af_lb  = CompactAFReader.initialize_from_file(filename_lb)
    cumulated_occupation_lb = af_lb.cumulated_occupation_array()
    occupation_lb = af_lb.occupation_array()


    af_hb  = CompactAFReader.initialize_from_file(filename_hb)
    cumulated_occupation_hb = af_hb.cumulated_occupation_array()
    occupation_hb = af_hb.occupation_array()
    #



    print("Coherent fraction EBS: ",cumulated_occupation_ebs[0])
    print("Coherent fraction LB: ",cumulated_occupation_lb[0])
    print("Coherent fraction HB: ",cumulated_occupation_hb[0])


    extensions = ["ebs","lb","hb"]
    data = [cumulated_occupation_ebs,cumulated_occupation_lb,cumulated_occupation_hb]
    data_occ = [occupation_ebs,occupation_lb,occupation_hb]


    plot(numpy.arange(cumulated_occupation_ebs.size),cumulated_occupation_ebs,
         numpy.arange(cumulated_occupation_lb.size),cumulated_occupation_lb,
         numpy.arange(cumulated_occupation_hb.size),cumulated_occupation_hb,
         legend=extensions)


    for i,extension in enumerate(extensions):
        f = open("cumulated_occupation_%s.dat"%extension,'w')
        data_i = data[i]
        for j in range(data_i.size):
            f.write("%d   %g  \n"%(j,data_i[j]))
        f.close()
        print("File written to disk: cumulated_occupation_%s.dat"%extension)

        f = open("occupation_%s.dat"%extension,'w')
        data_i = data_occ[i]
        for j in range(data_i.size):
            f.write("%d   %g  \n"%(j,data_i[j]))
        f.close()
        print("File written to disk: occupation_%s.dat"%extension)

    #
    # get indices
    #

    # first propagate a few modes only to check there are no errors
    # afp = AFpropagated.propagate(af,distance=distance,index_max=1,zoom=zoom)

