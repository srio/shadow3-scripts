import numpy
from dabam import dabam

def convert_daniele():
    input = "TXI_Flat1.dat"


    #
    # txt = "http://ftp.esrf.eu/pub/scisoft/dabam/data/dabam-067.dat"

    dm = dabam.initialize_from_external_data(input,
                               column_index_abscissas=0,
                               column_index_ordinates=3,
                               skiprows=1,
                               useHeightsOrSlopes=0,
                               to_SI_abscissas=1e-3,
                               to_SI_ordinates=1e-9,
                               detrending_flag=-1)

    if False:
        dm.metadata_set_info( YEAR_FABRICATION=None,
                              SURFACE_SHAPE=None,
                              FUNCTION=None,
                              LENGTH=None,
                              WIDTH=None,
                              THICK=None,
                              LENGTH_OPTICAL=None,
                              SUBSTRATE=None,
                              COATING=None,
                              FACILITY="LCLS",
                              INSTRUMENT=None,
                              POLISHING=None,
                              ENVIRONMENT=None,
                              SCAN_DATE=None,
                              CALC_HEIGHT_RMS=None,
                              CALC_HEIGHT_RMS_FACTOR=None,
                              CALC_SLOPE_RMS=None,
                              CALC_SLOPE_RMS_FACTOR=None,
                              USER_EXAMPLE=None,
                              USER_REFERENCE=input,
                              USER_ADDED_BY="Daniele Cocco and Manuel Sanchez del Rio",)

    dm.plot("heights")


    dm.write_output_dabam_files(filename_root="DABAM-%03d"%1,loaded_from_file=input)


def convert_antoine():
    from srxraylib.plot.gol import plot

    directory = "/home/manuel/google-drive/optics_COSMIC/"
    file = "COSMIC M101 LTP slope and height.csv"

    FILES = ["COSMIC_M101_LTP_slope_and_height.csv",
             "COSMIC_M111_LTP_slope_and_height.csv",
             "COSMIC_M112_LTP_slope_and_height.csv",
             "COSMIC_M121_LTP_slope_and_height.csv",
             "qerlinM201.dat",
             "QERLIN_M203_data_Trace_Line_1.csv",
             "amberM202faceUpResidSlope.dat"]

    for i,ifile in enumerate(FILES):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> file",ifile,type(file))
        # "amberM202faceUpResidSlope.dat"
        input = directory + ifile


        f = open(input,"r")
        txt = f.read()
        f.close()

        txt = txt.replace("g", " ")
        txt = txt.replace(",,,",",")
        txt = txt.replace(",,,",",")
        txt = txt.replace('"','')
        txt = txt.replace(",", "   ")

        txt = txt.split("\n")
        for ii in range(5):
            print(">>>>",txt[ii])


        a = numpy.loadtxt(txt, skiprows=1)
        print(a.shape)
        # plot(a[:,0],a[:,-1]*1e3,title="raw heights "+file)

        print(">> RAW slope error: %f urad RMS"%(a[:,1].std()*1e6))

        dm = dabam()

        dm = dabam.initialize_from_external_data(txt,
                                   column_index_abscissas=0,
                                   column_index_ordinates=1,
                                   skiprows=1,
                                   useHeightsOrSlopes=1,
                                   to_SI_abscissas=1e-3,
                                   to_SI_ordinates=1.0,
                                   detrending_flag=1)

        if True:

            if i <= 3:
                CALC_HEIGHT_RMS = a[:, 5].std() * 1e6
                CALC_HEIGHT_RMS_FACTOR = 1e-9
                CALC_SLOPE_RMS = a[:, 3].std() * 1e6
                CALC_SLOPE_RMS_FACTOR = 1e-6
                print(">> RAW slope error (detrended): %f urad RMS" % (a[:, 3].std() * 1e6))
                print(">> RAW height error (detrended): %f nm RMS" % (a[:, 5].std() * 1e3))
            else:
                CALC_HEIGHT_RMS        = None
                CALC_HEIGHT_RMS_FACTOR = None
                CALC_SLOPE_RMS         = None
                CALC_SLOPE_RMS_FACTOR  = None

            dm.metadata_set_info( YEAR_FABRICATION=None,
                                  SURFACE_SHAPE=None,
                                  FUNCTION=None,
                                  LENGTH=None,
                                  WIDTH=None,
                                  THICK=None,
                                  LENGTH_OPTICAL=(a[:,0].max()-a[:,0].min())*1e-3,
                                  SUBSTRATE=None,
                                  COATING=None,
                                  FACILITY="ALS",
                                  INSTRUMENT=None,
                                  POLISHING=None,
                                  ENVIRONMENT=None,
                                  SCAN_DATE=None,
                                  CALC_HEIGHT_RMS        = CALC_HEIGHT_RMS       ,
                                  CALC_HEIGHT_RMS_FACTOR = CALC_HEIGHT_RMS_FACTOR,
                                  CALC_SLOPE_RMS         = CALC_SLOPE_RMS        ,
                                  CALC_SLOPE_RMS_FACTOR  = CALC_SLOPE_RMS_FACTOR ,
                                  USER_EXAMPLE=None,
                                  USER_REFERENCE="From file: %s"%file,
                                  USER_ADDED_BY="Antoine Wojdyla and Manuel Sanchez del Rio",)

        # dm.plot("heights")
        #
        #

        print(dm.metadata)
        dm.write_output_dabam_files(filename_root="dabam-%03d"%(i+1),loaded_from_file=txt)

        print(">> RAW slope error: %f urad RMS"%(a[:,1].std()*1e6))
        if i <= 3:
            print(">> RAW slope error (detrended): %f urad RMS" % (a[:, 3].std() * 1e6))
            print(">> RAW height error (detrended): %f nm RMS" % (a[:, 5].std() * 1e6))

        print(">> DABAM slope error: %f urad RMS"%(dm.zSlopesUndetrended.std()*1e6))
        print(">> DABAM slope error (detrended): %f urad RMS"%(dm.zSlopes.std()*1e6))
        print(">> RAW height error (detrended): %f nm RMS" % (dm.zHeights.std() * 1e9))

    return i
if __name__ == "__main__":

    imax = convert_antoine()

    print("imax: ",imax)


