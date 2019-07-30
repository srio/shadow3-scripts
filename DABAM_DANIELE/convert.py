import numpy
from dabam import dabam
from srxraylib.plot.gol import plot
import matplotlib.pylab as plt

#


def load_external_profile(filename,
                     column_index_abscisa,
                     column_index_ordinates,
                     skiprows=1,
                     column_description=None,
                          type=None,  # "heights" or "slopes" -- deprecated, use FILE_FORMAT instead
                          FILE_FORMAT=2,  # 1=slopes, 2=heights
                          FILE_HEADER_LINES=0,
                          X1_FACTOR=1.0,
                          Y1_FACTOR=1.0,

                          YEAR_FABRICATION=None,
                          SURFACE_SHAPE="",
                          FUNCTION=None,
                          LENGTH=None,
                          WIDTH=None,
                          THICK=None,
                          LENGTH_OPTICAL=None,
                          SUBSTRATE=None,
                          COATING=None,
                          FACILITY=None,
                          INSTRUMENT=None,
                          POLISHING=None,
                          ENVIRONMENT=None,
                          SCAN_DATE=None,
                          # PLOT_TITLE_X1=None,
                          # PLOT_TITLE_Y1=None,
                          # PLOT_TITLE_Y2=None,
                          # PLOT_TITLE_Y3=None,
                          # PLOT_TITLE_Y4=None,
                          CALC_HEIGHT_RMS=None,
                          CALC_HEIGHT_RMS_FACTOR=None,
                          CALC_SLOPE_RMS=None,
                          CALC_SLOPE_RMS_FACTOR=None,
                          USER_EXAMPLE=None,
                          USER_REFERENCE=None,
                          USER_ADDED_BY=None,
                          file_out_root=None,

                          ):
    """
    load a profile from python arrays
    :param x: the abscissas (usually in m)
    :param y: the ordinates (heights in m or slopes in rad)
    :param type: set to 'heights' (default) of 'slopes' depending the type of input
    :return:
    """

    dm = dabam()
    a = numpy.loadtxt(filename,skiprows=skiprows)
    dm.rawdata = a

    dm.set_input_useAbscissasColumn(column_index_abscisa)
    dm.set_input_useOrdinatesColumn(column_index_ordinates)


    dm.set_input_entryNumber(-1)
    dm.set_input_multiply(1.0)
    dm.set_input_oversample(0.0)
    dm.set_input_setDetrending(-1)
    if FILE_FORMAT == 1:
        dm.set_input_useHeightsOrSlopes(1)
    elif FILE_FORMAT == 2:
        dm.set_input_useHeightsOrSlopes(0)

    dm.set_input_localFileRoot("<none (from python variable)>")



    dm.metadata = {}
    dm.metadata["FILE_FORMAT"] = FILE_FORMAT

    # this is redundant, now deprecated
    if type is not None:
        if type == "slopes":
            dm.metadata["FILE_FORMAT"] = 1
            dm.set_input_useHeightsOrSlopes(1)
        elif type == "heights":
            dm.metadata["FILE_FORMAT"] = 2
            dm.set_input_useHeightsOrSlopes(0)

    dm.metadata["FILE_HEADER_LINES"] = FILE_HEADER_LINES
    dm.metadata["X1_FACTOR"] = X1_FACTOR

    for i in range(1, a.shape[1]):
        dm.metadata["Y%d_FACTOR"%i] = Y1_FACTOR

    dm.metadata["Y_COLUMN_INDEX"] = column_index_ordinates


    dm.metadata["YEAR_FABRICATION"] = YEAR_FABRICATION
    dm.metadata["SURFACE_SHAPE"] = SURFACE_SHAPE
    dm.metadata["FUNCTION"] = FUNCTION
    dm.metadata["LENGTH"] = LENGTH
    dm.metadata["WIDTH"] = WIDTH
    dm.metadata["THICK"] = THICK
    dm.metadata["LENGTH_OPTICAL"] = LENGTH_OPTICAL
    dm.metadata["SUBSTRATE"] = SUBSTRATE
    dm.metadata["COATING"] = COATING
    dm.metadata["FACILITY"] = FACILITY
    dm.metadata["INSTRUMENT"] = INSTRUMENT
    dm.metadata["POLISHING"] = POLISHING
    dm.metadata["ENVIRONMENT"] = ENVIRONMENT
    dm.metadata["SCAN_DATE"] = SCAN_DATE
    if column_description is None:
        dm.metadata["PLOT_TITLE_X1"] = None
        for i in range(1, a.shape[1]):
            dm.metadata["PLOT_TITLE_Y%d"%i] = None
    else:
        try:
            dm.metadata["PLOT_TITLE_X1"] = column_description[0]
            for i in range(1,a.shape[1]+1):
                dm.metadata["PLOT_TITLE_X%d"] = column_description[i]
        except:
            dm.metadata["PLOT_TITLE_X1"] = None
            for i in range(1, a.shape[1]+1):
                dm.metadata["PLOT_TITLE_Y%d"] = None

    dm.metadata["CALC_HEIGHT_RMS"] = CALC_HEIGHT_RMS
    dm.metadata["CALC_HEIGHT_RMS_FACTOR"] = CALC_HEIGHT_RMS_FACTOR
    dm.metadata["CALC_SLOPE_RMS"] = CALC_SLOPE_RMS
    dm.metadata["CALC_SLOPE_RMS_FACTOR"] = CALC_SLOPE_RMS_FACTOR
    dm.metadata["USER_EXAMPLE"] = USER_EXAMPLE
    dm.metadata["USER_REFERENCE"] = USER_REFERENCE
    dm.metadata["USER_ADDED_BY"] = USER_ADDED_BY




    # self.metadata["SURFACE_SHAPE"] = ""
    # self.metadata["FACILITY"] = ""
    # self.metadata["CALC_SLOPE_RMS"] = None
    # self.metadata["CALC_HEIGHT_RMS"] = None

    # plot(dm.rawdata[:, column_index_abscisa], dm.rawdata[:, column_index_ordinates])

    for key in dm.metadata.keys():
        print(">>>>>>>>>>>>>>>>>>>> key",key)

    dm._make_calculations()

    dm.plot("heights")

    return dm

if __name__ == "__main__":
    filename = "TXI_Flat1.dat"

    # a  = numpy.loadtxt(filename,skiprows=1)
    #
    # print(a.shape)
    #
    # plot(a[:,0],a[:,3])
    #
    # dm = dabam()
    # dm.set_input_useAbscissasColumn(0)
    # dm.set_input_useOrdinatesColumn(1)

    dm = convert_to_dabam(filename,0,1,skiprows=1,
                             type=None,  # "heights" or "slopes" -- deprecated, use FILE_FORMAT instead
                             FILE_FORMAT=2, # 1 slopes in Col2
                                            # 2 = heights in Col2
                                            # 3 = slopes in Col2, file X1 Y1 X2 Y2
                                            # 4 = heights in Col2, file X1 Y1 X2 Y2
                             FILE_HEADER_LINES=0,
                             X1_FACTOR=1.0e-3,
                             Y1_FACTOR=1.0e-9,
                             YEAR_FABRICATION=None,
                             SURFACE_SHAPE="",
                             FUNCTION=None,
                             LENGTH=None,
                             WIDTH=None,
                             THICK=None,
                             LENGTH_OPTICAL=None,
                             SUBSTRATE=None,
                             COATING=None,
                             FACILITY=None,
                             INSTRUMENT=None,
                             POLISHING=None,
                             ENVIRONMENT=None,
                             SCAN_DATE=None,
                             # PLOT_TITLE_X1=None,
                             # PLOT_TITLE_Y1=None,
                             # PLOT_TITLE_Y2=None,
                             # PLOT_TITLE_Y3=None,
                             # PLOT_TITLE_Y4=None,
                             CALC_HEIGHT_RMS=None,
                             CALC_HEIGHT_RMS_FACTOR=None,
                             CALC_SLOPE_RMS=None,
                             CALC_SLOPE_RMS_FACTOR=None,
                             USER_EXAMPLE=None,
                             USER_REFERENCE=None,
                             USER_ADDED_BY=None,
                             )


    # dm.inputs["setDetrending"] = -1

    # print(dm.info_profiles())

    dm._write_output_dabam_files()