

import numpy
from collections import OrderedDict

from srundplug import calc3d_pysru, compare_flux
from srxraylib.plot.gol import plot_image, plot_surface, plot, plot_table


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

def plot_from_file(beamline_name,path="/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/"):
    filename = path+beamline_name+".npy"
    e,h,v,i = numpy.load(filename)

    for ie in range(e.size):
        plot_image(i[ie],h,v)

    return e,h,v,i

def plot_H_from_file(beamline_name,path="/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/",show=True):
    filename = path+beamline_name+".npy"
    e,h,v,i = numpy.load(filename)

    print("Shapes: e,h,v,i",e.shape,h.shape,v.shape,i.shape)
    tmp = i.sum(axis=2)
    plot_image(tmp,e,h,aspect='auto',cmap='Greys_r',xtitle='Photon energy (eV)',ytitle='Horizontal position (mm)',title=beamline_name,show=show)

    return e,h,v,i

def plot_H1D_from_file(beamline_name,lines=10,normalize=False,path="/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/",show=True):
    filename = path+beamline_name+".npy"
    e,h,v,i = numpy.load(filename)

    print("Shapes: e,h,v,i",e.shape,h.shape,v.shape,i.shape)
    tmp = i.sum(axis=2)
    tmp2 = numpy.zeros((lines,h.size))
    legend = []
    for i in range(lines):
        energy_index = int( i * (e.size/lines) )
        energy_value = e[energy_index]
        legend.append("E=%d eV"%energy_value)
        tmp2[i,:] = tmp[energy_index]
        if normalize: tmp2[i,:] /= tmp2[i,:].sum()

    plot_table(h,tmp2,xtitle='Horizontal position (mm)',ytitle='Intensity',legend=legend,show=show)

    return e,h,v,i

def plot_E_from_file(beamline_names,path="/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/",show=True):
    for index,beamline_name in enumerate(beamline_names):
        filename = path+beamline_name+".npy"
        e,h,v,i = numpy.load(filename)

        if index == 0:
            tmp = numpy.zeros((e.size,len(beamline_names)))
        print("Shapes: e,h,v,i",e.shape,h.shape,v.shape,i.shape)
        tmp[:,index] = i.sum(axis=(1,2))
    plot_table(e,tmp.T,xtitle='Photon energy (eV)',ytitle='Intensity', legend=beamline_names, show=show)


def plot_flux_from_file(beamline_name,path="/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/"):
    filename = path+beamline_name+".npy"
    e,h,v,i = numpy.load(filename)

    tmp = i.sum(axis=2)
    plot(e,tmp.sum(axis=1),xtitle='Photon energy (eV)',ytitle='Intensity',title=beamline_name,show=True)

    return e,h,v,i

def create_hdf5(mypath="/tmp_14_days/srio/ID24_EBS/"):
    import os,numpy
    file_list = os.listdir(mypath)

    for file in file_list:
        print("reading file: %s"%file)
        # i(e,h,v) stack e,h,v are 1D arrays; e: energy in eV, h,v: horizontal and vertical positions in mm
        e,h,v,i = numpy.load(file)
        # calculate i2(e,h)
        i2 = i.sum(axis=2)





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

    # calculate_beamline_3d("ID24",             photonEnergyMin=5000,photonEnergyMax=7500,photonEnergyPoints=501,save=True,zero_emittance=False)
    # calculate_beamline_3d("ID24_EBS",         photonEnergyMin=5000,photonEnergyMax=7500,photonEnergyPoints=501,save=True,zero_emittance=False)
    # calculate_beamline_3d("ID24_NO_EMITTANCE",photonEnergyMin=5000,photonEnergyMax=7500,photonEnergyPoints=501,save=True,zero_emittance=True)
    # calculate_beamline_3d("ID24_U32",         photonEnergyMin=5000,photonEnergyMax=7500,photonEnergyPoints=501,save=True,zero_emittance=False)
    # calculate_beamline_3d("ID24_EBS_U32",     photonEnergyMin=5000,photonEnergyMax=7500,photonEnergyPoints=501,save=True,zero_emittance=False)

    # plot_from_file("SAKURA-RUN1/ID24")
    # plot_from_file("SAKURA-RUN1/ID24_EBS")
    # plot_from_file("SAKURA-RUN1/ID24_NO_EMITTANCE")
    # plot_from_file("SAKURA-RUN1/ID24_U32")
    # plot_from_file("SAKURA-RUN1/ID24_U32")
    # plot_from_file("SAKURA-RUN1/ID24_EBS_U32")

    # plot_flux_from_file("SAKURA-RUN1/ID24_EBS")

    # plot_E_from_file(["ID24_EBS","ID24","ID24_NO_EMITTANCE"],show=True)

    # plot_H1D_from_file("ID24_EBS",lines=1,normalize=False,show=True)

    plot_H_from_file("ID24_EBS",show=False)
    plot_H_from_file("ID24",show=False)
    plot_H_from_file("ID24_NO_EMITTANCE",show=True)



    #
    # compare_flux(get_beamline("ID24_EBS"),emin=5000,emax=7500,include_pysru=False,zero_emittance=False,iplot=True)
    # compare_flux(get_beamline("ID24"),emin=5000,emax=7500,include_pysru=False,zero_emittance=False,iplot=True)
    # compare_flux(get_beamline("ID24_NO_EMITTANCE"),emin=5000,emax=7500,include_pysru=False,zero_emittance=False,iplot=True)

    # create_hdf5()