import numpy
from srundplug import calc3d_pysru, calc3d_srw, calc3d_urgent, calc3d_us

# #Physical constants (global, by now)
import scipy.constants as codata
codata_mee = numpy.array(codata.physical_constants["electron mass energy equivalent in MeV"][0])
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)





########################################################################################################################
#
# TEST XOPPY DEFAULTS
#
########################################################################################################################
def test_xoppy_calc_undulator_radiation(ELECTRONENERGY=6.04,ELECTRONENERGYSPREAD=0.001,ELECTRONCURRENT=0.2,\
                                       ELECTRONBEAMSIZEH=0.000395,ELECTRONBEAMSIZEV=9.9e-06,\
                                       ELECTRONBEAMDIVERGENCEH=1.05e-05,ELECTRONBEAMDIVERGENCEV=3.9e-06,\
                                       PERIODID=0.018,NPERIODS=222,KV=1.68,DISTANCE=30.0,GAPH=0.003,GAPV=0.003,\
                                       HSLITPOINTS=41,VSLITPOINTS=41,METHOD=3):
    print("Inside xoppy_calc_undulator_radiation. ")

    bl = {}
    bl['name'] = "test_xoppy_undulator_radiation"
    bl['ElectronBeamDivergenceH'] = ELECTRONBEAMDIVERGENCEH
    bl['ElectronBeamDivergenceV'] = ELECTRONBEAMDIVERGENCEV
    bl['ElectronBeamSizeH'] = ELECTRONBEAMSIZEH
    bl['ElectronBeamSizeV'] = ELECTRONBEAMSIZEV
    bl['ElectronCurrent'] = ELECTRONCURRENT
    bl['ElectronEnergy'] = ELECTRONENERGY
    bl['ElectronEnergySpread'] = ELECTRONENERGYSPREAD
    bl['Kv'] = KV
    bl['NPeriods'] = NPERIODS
    bl['PeriodID'] = PERIODID
    bl['distance'] = DISTANCE
    bl['gapH'] = GAPH
    bl['gapV'] = GAPV


    gamma = ELECTRONENERGY / (codata_mee * 1e-3)
    print ("Gamma: %f \n"%(gamma))

    resonance_wavelength = (1 + bl['Kv']**2 / 2.0) / 2 / gamma**2 * bl["PeriodID"]
    m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)
    resonance_energy = m2ev / resonance_wavelength

    print ("Resonance wavelength [A]: %g \n"%(1e10*resonance_wavelength))
    print ("Resonance energy [eV]: %g \n"%(resonance_energy))

    energy = None
    if energy == None:
        energy = resonance_energy+1


    outFile = "undulator_radiation.spec"


    if METHOD == 0:
        print("Undulator radiation calculation using US. Please wait...")
        e,h,v,p = calc3d_us(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
                                    photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
        print("Done")
    if METHOD == 1:
        print("Undulator radiation calculation using URGENT. Please wait...")
        e,h,v,p = calc3d_urgent(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
                                    photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
        print("Done")
    if METHOD == 2:
        print("Undulator radiation calculation using SRW. Please wait...")
        e,h,v,p = calc3d_srw(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
                                    photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
        print("Done")
    if METHOD == 3:
        print("Undulator radiation calculation using SRW. Please wait...")
        e,h,v,p = calc3d_pysru(bl,fileName=outFile,fileAppend=False,hSlitPoints=HSLITPOINTS,vSlitPoints=VSLITPOINTS,
                                    photonEnergyMin=energy,photonEnergyMax=13000,photonEnergyPoints=1,zero_emittance=False)
        print("Done")


    from srxraylib.plot.gol import plot_image
    plot_image(p[0],h,v)



########################################################################################################################
#
# Main code
#
########################################################################################################################

if __name__ == '__main__':



    #
    # XOPPY tests
    #

    test_xoppy_calc_undulator_radiation()
