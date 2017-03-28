__author__    = "Manuel Sanchez del Rio"
__contact__   = "srio@esrf.eu"
__copyright__ = "ESRF, 2014-2017"


# try:
#     raise Exception("use local srundplug")
#     from orangecontrib.xoppy.util import srundplug
# except:
import srundplug


import numpy
import scipy.constants as codata





########################################################################################################################
#
# Main code
#
########################################################################################################################

if __name__ == '__main__':



    read_dictionary = numpy.load('ID21noEmittance.npy').item()
    srundplug.calculate_power(read_dictionary)
    # srundplug.compare_flux_plot(read_dictionary,show=True)
    srundplug.compare_power_density_plot(read_dictionary,show=True)
    # srundplug.compare_radiation_plot(read_dictionary,show=True)


