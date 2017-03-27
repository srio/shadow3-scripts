__author__    = "Manuel Sanchez del Rio"
__contact__   = "srio@esrf.eu"
__copyright__ = "ESRF, 2014-2017"


# try:
#     raise Exception("use local srundplug")
#     from orangecontrib.xoppy.util import srundplug
# except:
import srundplug
from srundplug_examples import get_beamline

import numpy
from collections import OrderedDict
import sys




########################################################################################################################
#
# GLOBAL NAMES
#
########################################################################################################################

# #Physical constants (global, by now)
import scipy.constants as codata
codata_mee = (codata.m_e * codata.c**2 / codata.e) * 1e-6 # electron mass energy equivalent in MeV
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)


########################################################################################################################
#
# Main code
#
########################################################################################################################

if __name__ == '__main__':

    # zero_emittance = False
    iplot = True



    #
    # open spec file
    #

    fileName = None

    if fileName is not None:
        scanCounter = 0
        f = open(fileName,"w")
        f.write("#F "+fileName+"\n")
        f.close()


    beamline_names = ["ID21"] # "XRAY_BOOKLET","ID16_NA","ESRF_NEW_OB","SHADOW_DEFAULT"]
    bl = get_beamline("ID21")
    print(bl)


    #
    # Info
    #

    # for beamline_name in beamline_names:
    #     print(beamline_info(get_beamline(beamline_name,zero_emittance=zero_emittance),distance=26.0))



    #
    # Power density
    #
    # compare_power_density(get_beamline("ID21"    ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("SHADOW_DEFAULT" ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("XRAY_BOOKLET"   ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("ID16_NA"        ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("EBS_OB"),zero_emittance=zero_emittance,iplot=iplot)

    bl = srundplug.compare_power_density(bl,zero_emittance=True,post_convolution=True)
    if iplot: srundplug.compare_power_density_plot(bl,contour=True,surface=False)

    srundplug.calculate_power(bl)


    # #
    # # dump file
    # #
    # numpy.save(bl["name"],bl)
    #
    # read_dictionary = numpy.load('ID21.npy').item()
    # print(read_dictionary.keys())
    # calculate_power(read_dictionary)
    # # srundplug.compare_radiation_plot(read_dictionary,show=True)


