import numpy
import Shadow
from srxraylib.sources import srfunc
import matplotlib
from matplotlib import pylab as plt

matplotlib.rcParams.update({'font.size': 8})

#Physical constants (global, by now)
try:
    import scipy.constants.codata
    codata = scipy.constants.codata.physical_constants

    codata_c, tmp1, tmp2 = codata["speed of light in vacuum"]
    codata_c = numpy.array(codata_c)

    codata_mee, tmp1, tmp2 = codata["electron mass energy equivalent in MeV"]
    codata_mee = numpy.array(codata_mee)

    codata_me, tmp1, tmp2 = codata["electron mass"]
    codata_me = numpy.array(codata_me)

    codata_h, tmp1, tmp2 = codata["Planck constant"]
    codata_h = numpy.array(codata_h)

    codata_ec, tmp1, tmp2 = codata["elementary charge"]
    codata_ec = numpy.array(codata_ec)
except ImportError:
    print("Failed to import scipy. Finding alternative ways.")
    codata_c = numpy.array(299792458.0)
    codata_mee = numpy.array(9.10938291e-31)
    codata_h = numpy.array(6.62606957e-34)
    codata_ec = numpy.array(1.602176565e-19)

m2ev = codata_c*codata_h/codata_ec      # lambda(m)  = m2eV / energy(eV)



def calc_wiggler_spectrum(ener_gev=6.0,e_min=100.0,e_max=100000.00,file_field="",output_file=""):

    (traj, pars) = srfunc.wiggler_trajectory(b_from=1,
                                             inData=file_field,
                                             nPer=1,
                                             nTrajPoints=501,
                                             ener_gev=ener_gev,
                                             per=None,
                                             kValue=None,
                                             trajFile="tmp.traj")

    x,y = srfunc.wiggler_spectrum(traj, enerMin=e_min, enerMax=e_max,nPoints=500, \
                     electronCurrent=0.2, outFile="", elliptical=False)

    #
    tmp = (numpy.vstack((x,y)))
    print(tmp.shape)
    numpy.savetxt(output_file,tmp.T)

    xx = numpy.array((5000.,10000,20000,40000,80000))
    return numpy.interp(xx,x,y)

def calc_bm_spectrum(e_gev=6.0,e_min=100.0,e_max=100000.00,output_file=""):
    # input for ESRF BM
    import scipy.constants.codata
    r_m = 3.3*e_gev/0.856     # magnetic radius in m
    i_a = 0.2                 # electron current in A

    # calculate critical energy in eV
    gamma = e_gev*1e3/codata_mee
    ec_m = 4.0*numpy.pi*r_m/3.0/numpy.power(gamma,3) # wavelength in m
    ec_ev = m2ev/ec_m

    print("Critical energy = %f eV"%ec_ev)
    print("Magnetic radius = %f m"%r_m)
    print("Gamma = %f "%gamma)
    print("mee = %f "%codata_mee)
    energy_ev = numpy.linspace(e_min,e_max,500) # photon energy grid
    f_psi = 0    # flag: full angular integration
    flux = srfunc.sync_ene(f_psi,energy_ev,ec_ev=ec_ev,polarization=0,  \
           e_gev=e_gev,i_a=i_a,hdiv_mrad=1.0, \
           psi_min=0.0, psi_max=0.0, psi_npoints=1)

    # for 2mrad
    flux *= 2

    tmp = (numpy.vstack((energy_ev,flux)))
    print(tmp.shape)
    numpy.savetxt(output_file,tmp.T)

    x = numpy.array((5000.,10000,20000,40000,80000))
    return numpy.interp(x,energy_ev,flux)



if __name__ == "__main__":

    devices = ["3P","2PA","2PB","1P"]

    f = open("table_flux.txt","w")

    for device in devices:
        file_field = "SW_"+device+"cut.txt"
        tmp = calc_wiggler_spectrum(file_field=file_field,output_file="spectrum-%s.txt"%device)
        f.write(" %s & %4.3g & %4.3g & %4.3g & %4.3g & %4.3g %s \n"%(device,tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],r"\\"))


    tmp = calc_bm_spectrum(output_file="spectrum-1BM.txt")
    f.write(" BM & %4.3g & %4.3g & %4.3g & %4.3g & %4.3g %s \n"%(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],r"\\"))

    f.close()
    print("File written to disk: table_flux.txt")
