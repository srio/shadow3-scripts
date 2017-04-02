__author__    = "Manuel Sanchez del Rio"
__contact__   = "srio@esrf.eu"
__copyright__ = "ESRF, 2014-2017"


# try:
#     raise Exception("use local srundplug")
#     from orangecontrib.xoppy.util import srundplug
# except:
import srundplug


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


    elif nameBeamline == "ID16_NA":
        if silent == False:
            print("Setting inputs for ESRF-ID16_NA")

        ebeam['ElectronBeamDivergenceH'] = 10.3e-06
        ebeam['ElectronBeamDivergenceV'] = 1.2e-06
        ebeam['ElectronBeamSizeH'] = 0.0004131
        ebeam['ElectronBeamSizeV'] = 3.4e-06
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.04
        idv['Kv'] = 4.0 # 0.82
        idv['NPeriods'] = 77
        idv['PeriodID'] = 0.026
        drift['distance'] = 20.0
        slit['gapH'] = 15*0.001 #0.001
        slit['gapV'] = 15*0.001 #0.001

    elif nameBeamline == "ESRF_LB":
        if silent == False:
            print("Setting inputs for ESRF Low Beta")

        ebeam['ElectronBeamDivergenceH'] = 88.3e-6
        ebeam['ElectronBeamDivergenceV'] = 3.8e-6
        ebeam['ElectronBeamSizeH'] = 57e-6
        ebeam['ElectronBeamSizeV'] = 10.3e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.04
        idv['Kv'] = 1.68
        idv['NPeriods'] = int(4.0/0.018)
        idv['PeriodID'] = 0.018
        drift['distance'] = 30.0
        slit['gapH'] = 0.03 # 0.001
        slit['gapV'] = 0.005 # 0.001

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


    elif nameBeamline == "ESRF_HB":
        if silent == False:
            print("Setting inputs for ESRF High Beta")

        ebeam['ElectronBeamDivergenceH'] = 10.5e-6
        ebeam['ElectronBeamDivergenceV'] = 3.9e-6
        ebeam['ElectronBeamSizeH'] = 395e-6
        ebeam['ElectronBeamSizeV'] = 9.9e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.04
        idv['Kv'] = 1.68
        idv['NPeriods'] = int(4.0/0.018)
        idv['PeriodID'] = 0.018
        drift['distance'] = 30.0
        slit['gapH'] = 0.001 #0.001
        slit['gapV'] = 0.001 #0.001

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



    elif nameBeamline == "EBS_OB":
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

    elif nameBeamline == "ID24":
        if silent == False:
            print("Setting inputs for ESRF ID24/current")

        ebeam['ElectronBeamDivergenceH'] = 10.3e-6
        ebeam['ElectronBeamDivergenceV'] = 1.2e-6
        ebeam['ElectronBeamSizeH'] = 387.8e-6
        ebeam['ElectronBeamSizeV'] = 3.5e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.0
        idv['Kv'] = 1.272
        idv['PeriodID'] = 0.027
        idv['NPeriods'] = int(3.2/idv['PeriodID'])
        drift['distance'] = 27.0
        slit['gapH'] = 0.0025
        slit['gapV'] = 0.0025

    elif nameBeamline == "ID24_EBS":
        if silent == False:
            print("Setting inputs for ESRF ID24/EBS")

        ebeam['ElectronBeamDivergenceH'] = 5.2e-6
        ebeam['ElectronBeamDivergenceV'] = 1.4e-6
        ebeam['ElectronBeamSizeH'] = 27.2e-6
        ebeam['ElectronBeamSizeV'] = 3.4e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.0
        idv['Kv'] = 1.272
        idv['PeriodID'] = 0.027
        idv['NPeriods'] = int(3.2/idv['PeriodID'])
        drift['distance'] = 27.0
        slit['gapH'] = 0.0025
        slit['gapV'] = 0.0025

    elif nameBeamline == "ID24_NO_EMITTANCE":
        if silent == False:
            print("Setting inputs for ESRF ID24/EBS")

        ebeam['ElectronBeamDivergenceH'] = 5.2e-6
        ebeam['ElectronBeamDivergenceV'] = 1.4e-6
        ebeam['ElectronBeamSizeH'] = 27.2e-6
        ebeam['ElectronBeamSizeV'] = 3.4e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.0
        idv['Kv'] = 1.272
        idv['PeriodID'] = 0.027
        idv['NPeriods'] = int(3.2/idv['PeriodID'])
        drift['distance'] = 27.0
        slit['gapH'] = 0.0025 #0.001
        slit['gapV'] = 0.0025 #0.001

    elif nameBeamline == "ID24_U32":
        if silent == False:
            print("Setting inputs for ESRF ID24/current")

        ebeam['ElectronBeamDivergenceH'] = 10.3e-6
        ebeam['ElectronBeamDivergenceV'] = 1.2e-6
        ebeam['ElectronBeamSizeH'] = 387.8e-6
        ebeam['ElectronBeamSizeV'] = 3.5e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.0
        idv['Kv'] = 1.026#272
        idv['PeriodID'] = 0.032
        idv['NPeriods'] = int(3.2/idv['PeriodID'])
        drift['distance'] = 27.0
        slit['gapH'] = 0.0025 #0.001
        slit['gapV'] = 0.0025 #0.001

    elif nameBeamline == "ID24_EBS_U32":
        if silent == False:
            print("Setting inputs for ESRF ID24/EBS")

        ebeam['ElectronBeamDivergenceH'] = 5.2e-6
        ebeam['ElectronBeamDivergenceV'] = 1.4e-6
        ebeam['ElectronBeamSizeH'] = 27.2e-6
        ebeam['ElectronBeamSizeV'] = 3.4e-6
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.0
        idv['Kv'] = 1.026#272
        idv['PeriodID'] = 0.032
        idv['NPeriods'] = int(3.2/idv['PeriodID'])
        drift['distance'] = 27.0
        slit['gapH'] = 0.0025 #0.001
        slit['gapV'] = 0.0025 #0.001

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

    elif nameBeamline == "ID21":
        if silent == False:
            print("Setting inputs for ESRF ID21/current")

        ebeam['ElectronBeamDivergenceH'] = 0.000107
        ebeam['ElectronBeamDivergenceV'] = 1.16e-06
        ebeam['ElectronBeamSizeH'] = 4.99e-05
        ebeam['ElectronBeamSizeV'] = 3.45e-06
        ebeam['ElectronEnergySpread'] = 0.001
        ebeam['ElectronCurrent'] = 0.2
        ebeam['ElectronEnergy'] = 6.037
        idv['Kv'] = 2.14
        idv['PeriodID'] = 0.042
        idv['NPeriods'] = 37
        drift['distance'] = 26.0
        slit['gapH'] = 0.002
        slit['gapV'] = 0.002

    else:
        raise Exception("This name (%s) does not correspond at any name for a beamline"%nameBeamline)







    # build the beamline as a merged dictionary
    # TODO: merge elements is dangerous (possible key duplication) and
    #       does not allow multiple identical elements. Find a better solution...
    bl = OrderedDict()
    bl.update({'name':nameBeamline})
    bl.update(ebeam)
    bl.update(idv)
    bl.update(drift)
    bl.update(slit)

    if zero_emittance:
        bl = set_beamline_zero_emittance(bl)

    #if silent == False:
    #    print ("\n\n-----------------------------------------------------")
    #    for i,j in bl.items():
    #        print ("%s = %s" % (i,j) )
    #    print ("-----------------------------------------------------\n\n")
    return bl

def set_beamline_zero_emittance(bl):

    bl['ElectronBeamDivergenceH'] = 1e-30
    bl['ElectronBeamDivergenceV'] = 1e-30
    bl['ElectronBeamSizeH']       = 1e-30
    bl['ElectronBeamSizeV']       = 1e-30
    bl['ElectronEnergySpread']    = 1e-30

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
            # values that depen on screen distance
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



#
#
# ########################################################################################################################
# #
# # TEST XOPPY DEFAULTS
# #
# ########################################################################################################################
# def test_xoppy_calc_undulator_radiation(ELECTRONENERGY=6.04,ELECTRONENERGYSPREAD=0.001,ELECTRONCURRENT=0.2,\
#                                        ELECTRONBEAMSIZEH=0.000395,ELECTRONBEAMSIZEV=9.9e-06,\
#                                        ELECTRONBEAMDIVERGENCEH=1.05e-05,ELECTRONBEAMDIVERGENCEV=3.9e-06,\
#                                        PERIODID=0.018,NPERIODS=222,KV=1.68,DISTANCE=30.0,GAPH=0.003,GAPV=0.003,\
#                                        HSLITPOINTS=41,VSLITPOINTS=41,METHOD=3):
#     print("Inside xoppy_calc_undulator_radiation. ")
#
#     bl = OrderedDict()
#     bl['ElectronBeamDivergenceH'] = ELECTRONBEAMDIVERGENCEH
#     bl['ElectronBeamDivergenceV'] = ELECTRONBEAMDIVERGENCEV
#     bl['ElectronBeamSizeH'] = ELECTRONBEAMSIZEH
#     bl['ElectronBeamSizeV'] = ELECTRONBEAMSIZEV
#     bl['ElectronCurrent'] = ELECTRONCURRENT
#     bl['ElectronEnergy'] = ELECTRONENERGY
#     bl['ElectronEnergySpread'] = ELECTRONENERGYSPREAD
#     bl['Kv'] = KV
#     bl['NPeriods'] = NPERIODS
#     bl['PeriodID'] = PERIODID
#     bl['distance'] = DISTANCE
#     bl['gapH'] = GAPH
#     bl['gapV'] = GAPV
#
#
#     gamma = ELECTRONENERGY / (codata_mee * 1e-3)
#     print ("Gamma: %f \n"%(gamma))
#
#     resonance_wavelength = (1 + bl['Kv']**2 / 2.0) / 2 / gamma**2 * bl["PeriodID"]
#     m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)
#     resonance_energy = m2ev / resonance_wavelength
#
#     print ("Resonance wavelength [A]: %g \n"%(1e10*resonance_wavelength))
#     print ("Resonance energy [eV]: %g \n"%(resonance_energy))
#
#     energy = None
#     if energy == None:
#         energy = resonance_energy+1
#
#
#     outFile = "undulator_radiation.spec"
#
#
#     if METHOD == 0:
#         print("Undulator radiation calculation using US. Please wait...")
#         e,h,v,p = calc3d_us(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
#                                     photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
#         print("Done")
#     if METHOD == 1:
#         print("Undulator radiation calculation using URGENT. Please wait...")
#         e,h,v,p = calc3d_urgent(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
#                                     photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
#         print("Done")
#     if METHOD == 2:
#         print("Undulator radiation calculation using SRW. Please wait...")
#         e,h,v,p = calc3d_srw(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
#                                     photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
#         print("Done")
#     if METHOD == 3:
#         print("Undulator radiation calculation using SRW. Please wait...")
#         e,h,v,p = calc3d_pysru(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
#                                     photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
#         print("Done")
#
#
#     from srxraylib.plot.gol import plot_image
#     plot_image(p[0],h,v)
#
# ########################################################################################################################
# #
# # Calculate 3d maps (ID24, etc)
# #
# ########################################################################################################################
# def calculate_beamline_3d(beamline_name,photonEnergyMin=6500,photonEnergyMax=7500,photonEnergyPoints=5,
#                           save=True,zero_emittance=False):
#     hSlitPoints = 101
#     vSlitPoints = 101
#     e,h,v,i = calc3d_pysru(get_beamline(beamline_name,zero_emittance=zero_emittance) ,zero_emittance=zero_emittance,
#                 photonEnergyMin=photonEnergyMin,photonEnergyMax=photonEnergyMax,photonEnergyPoints=photonEnergyPoints,
#                 hSlitPoints=hSlitPoints,vSlitPoints=vSlitPoints,
#                 fileName="%s.spec"%beamline_name,fileAppend=False)
#     if save:
#         numpy.save("%s"%beamline_name,[e,h,v,i])
#
# def plot_from_file(beamline_name):
#     e,h,v,i = numpy.load("%s.npy"%beamline_name)
#     from srxraylib.plot.gol import plot_image
#     for ie in range(e.size):
#         plot_image(i[ie],h,v)
#
# def plot_H_from_file(beamline_name):
#     e,h,v,i = numpy.load("%s.npy"%beamline_name)
#     from srxraylib.plot.gol import plot_image, plot_surface, plot
#     tmp = i.sum(axis=2)
#     print(">>>",i.shape,tmp.shape)
#     plot_image(tmp,e,h,aspect='auto',cmap='Greys_r',xtitle='Photon energy (eV)',ytitle='Horizontal position (mm)',title=beamline_name,show=False)
#     plot(e,tmp.sum(axis=1),xtitle='Photon energy (eV)',ytitle='Intensity',title=beamline_name)
#


########################################################################################################################
#
# Main codes
#
########################################################################################################################

def k_scan():

    from srxraylib.plot.gol import plot,plot_show

    bl = get_beamline("EBS_OB")
    srundplug.USE_URGENT = False
    srundplug.USE_US = False


    print(bl)

    K_values = 1.86 * numpy.array([0.5,1.0,3.0])

    for K in K_values:
        B0 = 1.0
        bl['Kv'] = K
        bl['PeriodID'] = bl['Kv'] / (93.36 * B0)
        bl['NPeriods'] = int(1.0/bl['PeriodID'])

        # bl['gapH'] = 0.01
        # bl['gapV'] = 0.01

        fileName="K_scan_%4.2f"%bl['Kv']

        e_s,f_s = srundplug.calc1d_srw(bl,photonEnergyMin=2000.0,photonEnergyMax=50000.0,
              photonEnergyPoints=10000,zero_emittance=False,fileName=fileName+".spec",fileAppend=False)

        print(">>>>>> K: %f, Period: %f, N:%d "%(bl['Kv'],bl['PeriodID'],bl['NPeriods']))
        print("Power from integral of SRW spectrum: %f W"%(f_s.sum()*1e3*codata.e*(e_s[1]-e_s[0])))
        bl["calc1d_srw"] = {"energy":e_s,"flux":f_s}
        numpy.save(fileName,bl)

    # make plot to png

    data = []
    legend = []
    key = "calc1d_srw"

    for K in K_values:
        beamline_dict = numpy.load("K_scan_%4.2f.npy"%(K)).item()

        if key in beamline_dict.keys():
            data.append(beamline_dict[key]["energy"])
            data.append(beamline_dict[key]["flux"])
            legend.append("K=%4.2f"%(K))


    fig = plot(data,title=beamline_dict['name'],show=False,legend=legend,ylog=True,
             xtitle="Photon energy [eV]",ytitle="Photons/s/0.1%bw")

    import matplotlib.pylab as plt
    plt.savefig('k_scan.png',bbox_inches='tight')

    # plot_show()


def main():


    zero_emittance = True
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


    #
    # Info
    #

    # for beamline_name in beamline_names:
    #     print(beamline_info(get_beamline(beamline_name,zero_emittance=zero_emittance),distance=26.0))


    #
    # Flux
    #


    srundplug.USE_SRWLIB = True

    # compare_flux(get_beamline("XRAY_BOOKLET"  ),emin=100, emax=900,  npoints=200,zero_emittance=True,iplot=iplot)

    # this is quite slow
    #  compare_flux(get_beamline("ID16_NA"       ),emin=3000,emax=20000,npoints=100,zero_emittance=zero_emittance,
    #              iplot=iplot,srw_max_harmonic_number=61)

    # compare_flux(get_beamline("ESRF_NEW_OB"   ),emin=6500,emax=9500, npoints=200,zero_emittance=zero_emittance,iplot=iplot)
    # compare_flux(get_beamline("SHADOW_DEFAULT"),emin=3000,emax=50000,npoints=200,zero_emittance=zero_emittance,iplot=iplot)
    # compare_flux_from_3d("ESRF_NEW_OB",emin=6500,emax=9500,npoints=200,zero_emittance=zero_emittance,iplot=iplot)


    bl = srundplug.compare_flux(bl,emin=1000, emax=100000,  npoints=1000,zero_emittance=zero_emittance)
    if iplot: srundplug.plot_flux(bl)
    print("K: %f, Period: %f, N:%d "%(bl['Kv'],bl['PeriodID'],bl['NPeriods']))


    #
    # Power density
    #
    # compare_power_density(get_beamline("ID21"    ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("SHADOW_DEFAULT" ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("XRAY_BOOKLET"   ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("ID16_NA"        ),zero_emittance=zero_emittance,iplot=iplot)
    # compare_power_density(get_beamline("EBS_OB"),zero_emittance=zero_emittance,iplot=iplot)

    # bl = srundplug.compare_power_density(bl,zero_emittance=zero_emittance)
    # if iplot: srundplug.compare_power_density_plot(bl)

    #
    # Radiance
    #

    # for beamline_name in beamline_names:
    #     compare_radiation(get_beamline(beamline_name,zero_emittance=zero_emittance),     energy=None,zero_emittance=zero_emittance,iplot=True,show=True)
    # srundplug.USE_PYSRU = True
    # bl = srundplug.compare_radiation(bl,
    #                                  photonEnergyMin=1000,photonEnergyMax=100000,photonEnergyPoints=500,
    #                                  zero_emittance=zero_emittance,iplot=True,show=True)


    if iplot: srundplug.plot_radiation(bl,show=True)






    #
    # dump file
    # #
    # numpy.save("ID21noEmittance",bl)

    # read_dictionary = numpy.load('ID21.npy').item()
    # print(read_dictionary.keys())
    # srundplug.calculate_power(read_dictionary)
    # srundplug.compare_power_density_plot(read_dictionary,show=True)
    # srundplug.compare_radiation_plot(read_dictionary,show=True)


def k_value(photon_energy, electron_energy=6.04, period=0.018, ):

    gamma = electron_energy* 1e9 / (codata.m_e *  codata.c**2 / codata.e)

    lambdan = codata.h*codata.c/codata.e*1e10 / photon_energy # in A

    harm_number = -1
    KK = -1
    while KK < 0:
        harm_number += 2
        lambda1 = lambdan * harm_number
        KK = ( 2*( (1e-10*lambda1)/(period) *2*gamma*gamma  - 1))

    return harm_number,numpy.sqrt(KK)

def tc():
    emin = 2000.0
    emax = 200000.0
    epoints = 100

    energies = numpy.linspace(emin,emax,epoints)

    for energy in energies:
        H,K = k_value(energy)
        print(">>>>E=%5.3f H:%d K:%5.3f "%(energy,H,K))



if __name__ == '__main__':
    # main()
    # k_scan()
    tc()