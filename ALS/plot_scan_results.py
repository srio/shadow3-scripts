
import numpy
from srxraylib.plot.gol import plot

def get_shadow_result(file_shadow_results,keyword):
    # file_shadow_results = "cosmic_scan_energy_grating806eV.dat"
    a = numpy.loadtxt(file_shadow_results,skiprows=4)

    # L photon energy [eV]  alpha [deg]  beta [deg]  size source V [um]  size [um]  source*M  C  resolving power

    energy = a[:,0]
    alpha = a[:,1]
    beta = a[:,2]
    source_size = a[:,3]
    image_size = a[:,4]
    source_size_times_magnification = a[:,5]
    footprint = a[:,6]
    C = a[:,7]
    resolving_power = a[:,8]

    if keyword == "alpha":
        return energy,alpha
    elif keyword == "beta":
        return energy,beta
    elif keyword == "source_size":
        return energy,source_size
    elif keyword == "image_size":
        return energy,image_size
    elif keyword == "source_size_times_magnification":
        return energy,source_size_times_magnification
    elif keyword == "footprint":
        return energy,footprint
    elif keyword == "C":
        return energy,C
    elif keyword == "resolving_power":
        return energy,resolving_power

def plot_shadow_results(filename,keyword):
    x,y = get_shadow_result(filename,keyword)
    plot(x,y,xtitle="Photon Energy [eV]",ytitle=keyword)
if __name__ == "__main__":
    filename = "cosmic_scan_energy_grating806eV.dat"
    plot_shadow_results(filename,"resolving_power")


