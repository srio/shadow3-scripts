import os
import numpy
from urllib.request import urlretrieve
from silx.io.specfile import SpecFile


def get_dabax_file(filename, url="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/"):

    try:
        if os.path.exists(filename):
            print("File exists: %s " % filename)
        else:
            filepath, http_msg = urlretrieve(url + filename,
                        filename=filename,
                        reporthook=None,
                        data=None)

            print("File %s downloaded from %s" % (filepath, url + filename))
        return True
    except:
        return False

#
# f0
#

def get_data_from_dabax_file(entry_name="Si", filename="MassEnergyAbsorption_NIST.dat"):

    # if getattr(get_f0_coeffs_from_dabax_file,'sf') is not None:
    #     sf = getattr(get_f0_coeffs_from_dabax_file,'sf')
    # else:
    error_flag = get_dabax_file(filename)
    if error_flag == False:
        raise(FileNotFoundError)
    sf = SpecFile(filename)
        # get_f0_coeffs_from_dabax_file.sf = sf

    flag_found = False

    for index in range(len(sf)):
        s1 = sf[index]
        name = s1.scan_header_dict["S"]

        if name.split(' ')[1] == entry_name:
            flag_found = True
            index_found = index

    if flag_found:
        return (sf[index_found].data)
    else:
        return []



if __name__ == "__main__":
    from srxraylib.plot.gol import plot, set_qt
    import xraylib
    set_qt()

    get_dabax_file("MassEnergyAbsorption_NIST.dat")

    for Z in [4, 6, 14, 29, 82]: # range(1, 93):
        symbol = xraylib.AtomicNumberToSymbol(Z)
        data = get_data_from_dabax_file(symbol)
        print(data.shape)
        energy = data[0, :]
        mu_en = data[2, :]

        mu_en_xraylib = numpy.zeros_like(mu_en)
        mu_xraylib = numpy.zeros_like(mu_en)

        for i, ienergy in enumerate(energy):
            mu_en_xraylib[i] = xraylib.CS_Energy(Z, 1e-3 * ienergy)
            try:
                mu_xraylib[i] = xraylib.CS_Total(Z, 1e-3 * ienergy)
            except:
                pass

        # plot(energy, mu_en,
        #      energy, mu_en_xraylib,
        #      xlog=1,ylog=1, title=symbol,legend=["MassEnergyAbsorption_NIST","CS_Energy xraylib"],
        #      xtitle="Photon energy [eV]", ytitle="Cross section cm2/g", xrange=[1e3,5e5], show=0)

        plot(1e-3 * energy, mu_xraylib,
             1e-3 * energy, mu_en,
             xlog=1,ylog=1, title=symbol,legend=["CS_Total xraylib","MassEnergyAbsorption_NIST"],
             xtitle="Photon energy [keV]", ytitle="Cross section cm2/g", xrange=[1,500], show=1)