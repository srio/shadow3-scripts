from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy

from srxraylib.util.h5_simple_writer import H5SimpleWriter
#
########################################################################################################################
#




if __name__ == "__main__":

    from srxraylib.plot.gol import plot
    import matplotlib.pylab as plt

    filename_ebs="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"


    af  = CompactAFReader.initialize_from_file(filename_ebs)

    # print(af.info())

    c = af.cumulated_occupation_array()


    f = plot(numpy.arange(c.size),c,xtitle="Mode index",ytitle="Cumulative occupation",yrange=[0,1],xlog=True,show=False)
    plt.savefig("/tmp/vx_cumulated.png")
    print("File written to disk: /tmp/vx_cumulated.png")

    plt.show()




