import numpy
from dabam import dabam


if __name__ == "__main__":

    input = "TXI_Flat1.dat"


    #
    # txt = "http://ftp.esrf.eu/pub/scisoft/dabam/data/dabam-067.dat"

    dm = dabam.initialize_from_external_data(input,
                               column_index_abscisa=0,
                               column_index_ordinates=3,
                               skiprows=1,
                               useHeightsOrSlopes=0,
                               to_SI_abscissas=1e-3,
                               to_SI_ordinates=1e-9,
                               detrending_flag=-1)

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


