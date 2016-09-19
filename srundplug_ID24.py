

import numpy
from collections import OrderedDict

from srundplug import calc3d_pysru
from srxraylib.plot.gol import plot_image, plot_surface, plot


__author__    = "Manuel Sanchez del Rio"
__contact__   = "srio@esrf.eu"
__copyright__ = "ESRF, 2016"

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

    if nameBeamline == "ID24":
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

    else:
        raise Exception("This name (%s) does not correspond at any name for a beamline"%nameBeamline)

    if zero_emittance:
        ebeam['ElectronBeamDivergenceH'] = 1e-30
        ebeam['ElectronBeamDivergenceV'] = 1e-30
        ebeam['ElectronBeamSizeH']       = 1e-30
        ebeam['ElectronBeamSizeV']       = 1e-30
        ebeam['ElectronEnergySpread']    = 1e-30

    bl = OrderedDict()
    bl.update({'name':nameBeamline})
    bl.update(ebeam)
    bl.update(idv)
    bl.update(drift)
    bl.update(slit)

    return bl


def calculate_beamline_3d(beamline_name,photonEnergyMin=6500,photonEnergyMax=7500,photonEnergyPoints=5,
                          save=False,zero_emittance=False):
    hSlitPoints = 101
    vSlitPoints = 101
    e,h,v,i = calc3d_pysru(get_beamline(beamline_name,zero_emittance=zero_emittance) ,zero_emittance=zero_emittance,
                photonEnergyMin=photonEnergyMin,photonEnergyMax=photonEnergyMax,photonEnergyPoints=photonEnergyPoints,
                hSlitPoints=hSlitPoints,vSlitPoints=vSlitPoints,
                fileName="%s.spec"%beamline_name,fileAppend=False)
    if save:
        numpy.save("%s"%beamline_name,[e,h,v,i])

def plot_from_file(beamline_name):
    e,h,v,i = numpy.load("%s.npy"%beamline_name)
    for ie in range(e.size):
        plot_image(i[ie],h,v)

def plot_H_from_file(beamline_name):
    e,h,v,i = numpy.load("%s.npy"%beamline_name)
    tmp = i.sum(axis=2)
    plot_image(tmp,e,h,aspect='auto',cmap='Greys_r',xtitle='Photon energy (eV)',ytitle='Horizontal position (mm)',title=beamline_name,show=True)

def plot_flux_from_file(beamline_name):
    e,h,v,i = numpy.load("%s.npy"%beamline_name)
    tmp = i.sum(axis=2)
    plot(e,tmp.sum(axis=1),xtitle='Photon energy (eV)',ytitle='Intensity',title=beamline_name,show=True)


########################################################################################################################
#
# Main code
#
########################################################################################################################

if __name__ == '__main__':

    zero_emittance = False
    iplot = True




    #
    # ID24
    #

    # calculate_beamline_3d("ID24",             photonEnergyMin=6500,photonEnergyMax=8000,photonEnergyPoints=501,save=True,zero_emittance=False)
    # calculate_beamline_3d("ID24_EBS",         photonEnergyMin=6500,photonEnergyMax=8000,photonEnergyPoints=501,save=True,zero_emittance=False)
    # calculate_beamline_3d("ID24_NO_EMITTANCE",photonEnergyMin=6500,photonEnergyMax=8000,photonEnergyPoints=501,save=True,zero_emittance=True)
    # calculate_beamline_3d("ID24_U32",         photonEnergyMin=6500,photonEnergyMax=8000,photonEnergyPoints=501,save=True,zero_emittance=False)
    # calculate_beamline_3d("ID24_EBS_U32",     photonEnergyMin=6500,photonEnergyMax=8000,photonEnergyPoints=501,save=True,zero_emittance=False)

    # plot_from_file("ID24")
    # plot_from_file("ID24_EBS")
    # plot_from_file("ID24_NO_EMITTANCE")
    # plot_from_file("ID24_U32")
    # plot_from_file("ID24_U32")
    # plot_from_file("ID24_EBS_U32")

    plot_flux_from_file("ID24_EBS")
    # plot_H_from_file("ID24_EBS")
    #
    # compare_flux("ID24_EBS",emin=3000,emax=15000,include_pysru=False,zero_emittance=False,iplot=True)
