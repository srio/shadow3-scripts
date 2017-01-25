

try:
    from orangecontrib.xoppy.util import srundplug
    from srundplug import compare_radiation, compare_flux, compare_power_density
except:
    from srundplug import compare_radiation, compare_flux, compare_power_density

__author__    = "Manuel Sanchez del Rio"
__contact__   = "srio@esrf.eu"
__copyright__ = "ESRF, 2014-2017"

import numpy
from collections import OrderedDict
import sys

#catch standard optput
try:
    from io import StringIO  # Python3
except ImportError:
    from StringIO import StringIO  # Python2

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
# Beamline parameters
#
########################################################################################################################
def get_beamline(nameBeamline,zero_emittance=False,silent=False):

    #
    # get the elements
    #

    ebeam = OrderedDict()
    idv   = OrderedDict()
    drift = OrderedDict()
    slit  = OrderedDict()

    #
    # modify elements
    #

    if nameBeamline == "XRAY_BOOKLET":
        if silent == False:
            print("Setting inputs for fig 2-5 in X-ray Data Booklet")

        ebeam['ElectronBeamDivergenceH'] = 1e-20
        ebeam['ElectronBeamDivergenceV'] = 1e-20
        ebeam['ElectronBeamSizeH'] = 1e-20
        ebeam['ElectronBeamSizeV'] = 1e-20
        ebeam['ElectronEnergySpread'] = 1e-20
        ebeam['ElectronCurrent'] = 1.0
        ebeam['ElectronEnergy'] = 1.3
        idv['Kv'] = 1.87
        idv['NPeriods'] = 14
        idv['PeriodID'] = 0.035
        drift['distance'] =   1.0*1e2
        slit['gapH']      = 0.002*1e2 #0.001
        slit['gapV']      = 0.002*1e2 #0.001



    elif nameBeamline == "ESRF_LB_OB":
        if silent == False:
            print("Setting inputs for ESRF Low Beta after OB")

        ebeam['ElectronBeamDivergenceH'] = 106.9e-6
        ebeam['ElectronBeamDivergenceV'] = 1.2e-6
        ebeam['ElectronBeamSizeH'] = 37.4e-6
        ebeam['ElectronBeamSizeV'] = 3.5e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.04
        idv['Kv'] = 1.68
        idv['NPeriods'] = int(4.0/0.018)
        idv['PeriodID'] = 0.018
        drift['distance'] = 30.0
        slit['gapH'] = 0.001
        slit['gapV'] = 0.001


    elif nameBeamline == "ESRF_HB_OB":
        if silent == False:
            print("Setting inputs for ESRF High Beta after OB")

        ebeam['ElectronBeamDivergenceH'] = 10.3e-6
        ebeam['ElectronBeamDivergenceV'] = 1.2e-6
        ebeam['ElectronBeamSizeH'] = 387.8e-6
        ebeam['ElectronBeamSizeV'] = 3.5e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.04
        idv['Kv'] = 1.68
        idv['NPeriods'] = int(4.0/0.018)
        idv['PeriodID'] = 0.018
        drift['distance'] = 30.0
        slit['gapH'] = 0.001
        slit['gapV'] = 0.001



    elif nameBeamline == "ESRF_NEW_OB":
        if silent == False:
            print("Setting inputs for ESRF New Lattice")

        ebeam['ElectronBeamDivergenceH'] = 5.2e-6
        ebeam['ElectronBeamDivergenceV'] = 1.4e-6
        ebeam['ElectronBeamSizeH'] = 27.2e-6
        ebeam['ElectronBeamSizeV'] = 3.4e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.0
        idv['Kv'] = 1.68
        idv['NPeriods'] = int(4.0/0.018)
        idv['PeriodID'] = 0.018
        drift['distance'] = 30.0
        slit['gapH'] = 0.001
        slit['gapV'] = 0.001

    elif nameBeamline == "SHADOW_DEFAULT":
        if silent == False:
            print("Setting inputs for SHADOW_DEFAULT")

        ebeam['ElectronBeamSizeH'] = 0.04e-2
        ebeam['ElectronBeamSizeV'] = 0.001e-2
        ebeam['ElectronBeamDivergenceH'] = 4e-9 / ebeam['ElectronBeamSizeH']
        ebeam['ElectronBeamDivergenceV'] = 4e-11 / ebeam['ElectronBeamSizeV']
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.04
        idv['Kv'] = 0.25
        idv['NPeriods'] = 50.0
        idv['PeriodID'] = 0.032
        drift['distance'] = 10.0
        slit['gapH'] = 2.0e-3
        slit['gapV'] = 2.0e-3


    elif nameBeamline == "ELETTRA":
        if silent == False:
            print("Setting inputs for Luca")
        # emittances from https://www.researchgate.net/profile/Emanuel_Karantzoulis/publication/277300230_EVOLUTION_OF_ELETTRA_TOWARDS_AN_ULTIMATE_LIGHT_SOURCE/links/5565e9ea08aeccd777359b37.pdf?origin=publication_list
        ebeam['ElectronBeamSizeH'] = 240e-6
        ebeam['ElectronBeamSizeV'] = 14e-6
        ebeam['ElectronBeamDivergenceH'] = 7e-9 / ebeam['ElectronBeamSizeH']
        ebeam['ElectronBeamDivergenceV'] = 0.01*7e-9 / ebeam['ElectronBeamSizeV']
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.4
        ebeam['ElectronEnergy'] = 2.0
        idv['Kv'] = 2.3
        idv['PeriodID'] = 0.025
        idv['NPeriods'] = int(1.5/idv['PeriodID'])
        drift['distance'] = 27.0
        slit['gapH'] = 0.0025 #0.001
        slit['gapV'] = 0.0025 #0.001
    else:
        raise Exception("This name (%s) does not correspond at any name for a beamline"%nameBeamline)

    if zero_emittance:
        ebeam['ElectronBeamDivergenceH'] = 1e-30
        ebeam['ElectronBeamDivergenceV'] = 1e-30
        ebeam['ElectronBeamSizeH']       = 1e-30
        ebeam['ElectronBeamSizeV']       = 1e-30
        ebeam['ElectronEnergySpread']    = 1e-30




    # build the beamline as a merged dictionary
    # TODO: merge elements is dangerous (possible key duplication) and
    #       does not allow multiple identical elements. Find a better solution...
    bl = OrderedDict()
    bl.update({'name':nameBeamline})
    bl.update(ebeam)
    bl.update(idv)
    bl.update(drift)
    bl.update(slit)

    return bl

def beamline_info(bl,photonEnergy=None,distance=None,silent=False):

    #init capture standard output
    # see http://wrongsideofmemphis.com/2010/03/01/store-standard-output-on-a-variable-in-python/
    old_stdout = sys.stdout
    result = StringIO()
    sys.stdout = result

    print("=============  Undulator parameters =============================\n")

    print("Inputs: \n")
    for i,j in bl.items():
        print ("%s = %s" % (i,j) )

    print("\n\nOutputs (supposing at waist): \n")
    print ("Electron beam Emittance H [m.rad]: %e \n"%(bl['ElectronBeamSizeH']*bl['ElectronBeamDivergenceH']))
    print ("Electron beam Emittance V [m.rad]: %e \n"%(bl['ElectronBeamSizeV']*bl['ElectronBeamDivergenceV']))
    if bl['ElectronBeamDivergenceH'] != 0.0:
        print ("Electron Beta H [m]: %e \n"%(bl['ElectronBeamSizeH']/bl['ElectronBeamDivergenceH']))
    if bl['ElectronBeamDivergenceV'] != 0.0:
        print ("Electron Beta V [m]: %e \n"%(bl['ElectronBeamSizeV']/bl['ElectronBeamDivergenceV']))
    l1 = bl['PeriodID']*bl['NPeriods']
    print ("Undulator length [m]: %f \n"%(l1))


    gamma = bl['ElectronEnergy'] / (codata_mee * 1e-3)
    print ("Gamma: %f \n"%(gamma))

    resonance_wavelength = (1 + bl['Kv']**2 / 2.0) / 2 / gamma**2 * bl["PeriodID"]
    resonance_energy = m2ev / resonance_wavelength

    print ("Resonance wavelength [A]: %g \n"%(1e10*resonance_wavelength))
    print ("Resonance energy [eV]: %g \n"%(resonance_energy))


    #
    # energy-dependent parameters
    #
    if photonEnergy is None:
        photonEnergy = [resonance_energy,3*resonance_energy,5*resonance_energy]

    if photonEnergy != None:
        phE = numpy.array(photonEnergy)
        phE.shape = -1
        for phEi in phE:
            print ("\n\n----------------------  photon energy [eV]: %0.2f \n"%(phEi))
            print('\n')
            lambda1 = m2ev/phEi
            print ("   photon wavelength [A]: %f \n"%(lambda1*1e10))

            # calculate sizes of the photon undulator beam
            # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
            s_phot = 2.740/(4e0*numpy.pi)*numpy.sqrt(l1*lambda1)
            sp_phot = 0.69*numpy.sqrt(lambda1/l1)
            print('\n')
            print('   RMS electon size H/V [um]: '+
                 repr(bl['ElectronBeamSizeH']*1e6)+ ' /  '+
                 repr(bl['ElectronBeamSizeV']*1e6) )
            print('   RMS electon divergence H/V[urad]: '+
                 repr(bl['ElectronBeamDivergenceH']*1e6)+ ' /  '+
                 repr(bl['ElectronBeamDivergenceV']*1e6)  )
            print('\n')
            print('   RMS radiation size [um]: '+repr(s_phot*1e6))
            print('   RMS radiation divergence [urad]: '+repr(sp_phot*1e6))
            print('\n')
            print('   Photon beam (convolution): ')
            photon_h = numpy.sqrt(numpy.power(bl['ElectronBeamSizeH'],2) + numpy.power(s_phot,2) )
            photon_v = numpy.sqrt(numpy.power(bl['ElectronBeamSizeV'],2) + numpy.power(s_phot,2) )
            photon_hp = numpy.sqrt(numpy.power(bl['ElectronBeamDivergenceH'],2) + numpy.power(sp_phot,2) )
            photon_vp = numpy.sqrt(numpy.power(bl['ElectronBeamDivergenceV'],2) + numpy.power(sp_phot,2) )
            print('   RMS size H/V [um]: '+ repr(photon_h*1e6) + '  /  '+repr(photon_v*1e6))
            print('   RMS divergence H/V [um]: '+ repr(photon_hp*1e6) + '  /  '+repr(photon_vp*1e6))

            print('\n')
            cohH = lambda1/4/numpy.pi / photon_h / photon_hp
            cohV = lambda1/4/numpy.pi / photon_v / photon_vp
            print('   Coherent volume in H phase space: '+ repr(cohH) )
            print('   Coherent volume in V phase space: '+ repr(cohV) )
            print('\n')
            dls = numpy.sqrt(2*l1*lambda1)/4/numpy.pi
            print('   RMS diffraction limit source size [um]: '+ repr(dls*1e6) )
            print('   FWHM diffraction limit source size [um]: '+ repr(dls*2.35*1e6) )
            #
            # values that depend on screen distance
            #
            if distance == None:
                distance = bl['distance']

            if distance != None:
                print('\n')
                hRMS = numpy.sqrt( numpy.power(photon_hp*distance,2) + numpy.power(photon_h,2))
                vRMS = numpy.sqrt( numpy.power(photon_vp*distance,2) + numpy.power(photon_v,2))

                print('   At a screen placed at :%f m from the source:\n'%(distance))
                print('   RMS size H/V [mm]: '+ repr(hRMS*1e3) + '  /  '+repr(vRMS*1e3))
                print('   FWHM size H/V [mm]: '+ repr(hRMS*2.35*1e3) + '  /  '+repr(vRMS*2.35*1e3))
                print('\n')
                print('   FWHM coherence length H [um] : '+ repr(hRMS*cohH*2.35*1e6) )
                print('   FWHM coherence length V [um] : '+ repr(vRMS*cohV*2.35*1e6) )


    print("=================================================================\n")

    sys.stdout = old_stdout
    result_string = result.getvalue()

    if silent == False:
        print(result_string)

    return  result_string


########################################################################################################################
#
# Main code
#
########################################################################################################################

if __name__ == '__main__':

    zero_emittance = False
    iplot = True
    include_pysru = False # otherwise flux and power_density are extremely slow

    #
    # open spec file to collect all results
    #

    fileName = None

    if fileName is not None:
        scanCounter = 0
        f = open(fileName,"w")
        f.write("#F "+fileName+"\n")
        f.close()

    beamline_names = ["ELETTRA"] # "XRAY_BOOKLET","ESRF_NEW_OB","SHADOW_DEFAULT"]

    #
    # Info
    #

    for beamline_name in beamline_names:
        print(beamline_info(get_beamline(beamline_name,zero_emittance=zero_emittance),distance=20.0))

    #
    # Radiance
    #

    for beamline_name in beamline_names:
        compare_radiation(get_beamline(beamline_name,zero_emittance=zero_emittance),     energy=None,zero_emittance=zero_emittance,iplot=True,show=True)


    #
    # Flux
    #

    compare_flux(get_beamline("ELETTRA"  ),emin=100, emax=2000,  npoints=200,zero_emittance=zero_emittance,include_pysru=include_pysru,iplot=iplot)


    #
    # Power density
    #

    compare_power_density(get_beamline("ELETTRA"),include_pysru=include_pysru,zero_emittance=zero_emittance,iplot=iplot)
