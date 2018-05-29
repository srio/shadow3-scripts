
import numpy
from orangecontrib.xoppy.util.xoppy_undulators import xoppy_calc_undulator_power_density, xoppy_calc_undulator_spectrum
from orangecontrib.xoppy.util.xoppy_xraylib_util import xpower_calc
from orangecontrib.xoppy.util.fit_gaussian2d import fit_gaussian2d, info_params, twoD_Gaussian

from srxraylib.plot.gol import plot, plot_image
import scipy.constants as codata

def calculate_line(photon_energy,undulator_period,N,K,thickness_diamond_mm,distance,slit_h,slit_v,coating,incident_angle_mrad,
                   do_plot=False):


    print("#########################   INPUTS  ###################################")
    print("photon_energy=",photon_energy)
    print("undulator_period=",undulator_period)
    print("N=",N)
    print("K=",K)
    print("thickness_diamond_mm=",thickness_diamond_mm)
    print("distance=",distance)
    print("slit_h=",slit_h)
    print("coating=",coating)
    print("incident_angle_mrad=",incident_angle_mrad)
    print("#######################################################################")

    out_dictionary = {}

    #
    # Spectrum simulation
    #

    #ULATTICEFILE S28D.mat
    #UEPSILONX      1.3166e-10
    #UEPSILONY           5e-12
    #BETAX    = 6.89997
    #BETAY    = 2.6447
    SIGMAX   = 30.1836 * 1e-6
    SIGMAY   = 3.63641 * 1e-6
    SIGMAXP  = 4.36821 * 1e-6
    SIGMAYP  = 1.37498 * 1e-6

    METHOD                   = 2       # US=0 URGENT=1 SRW=2


    print("\n\n Computing spectrum \n\n")
    e, f, spectral_power, cumulated_power = \
        xoppy_calc_undulator_spectrum(ELECTRONENERGY=6.0,ELECTRONENERGYSPREAD=0.001,ELECTRONCURRENT=0.2,\
                                  ELECTRONBEAMSIZEH=SIGMAX,ELECTRONBEAMSIZEV=SIGMAY,\
                                  ELECTRONBEAMDIVERGENCEH=SIGMAXP,ELECTRONBEAMDIVERGENCEV=SIGMAYP,\
                                  PERIODID=undulator_period,NPERIODS=N,KV=K,DISTANCE=distance,GAPH=slit_h,GAPV=slit_v,\
                                  PHOTONENERGYMIN=1000.0,PHOTONENERGYMAX=100000.0,PHOTONENERGYPOINTS=500,METHOD=2,
                                  USEEMITTANCES=1)


    power_in_spectrum = f.sum()*1e3*codata.e*(e[1]-e[0])

    out_dictionary["power_in_spectrum"] = power_in_spectrum
    if do_plot:
        plot(e,spectral_power,title="E = %d keV"%photon_energy)

    #
    # optical system
    #

    # """
    # Apply reflectivities/transmittivities of optical elements on a source spectrum
    #
    # :param energies: the array with photon energies in eV
    # :param source: the spectral intensity or spectral power
    # :param substance: a list with descriptors of each optical element  material
    # :param flags: a list with 0 (filter or attenuator) or 1 (mirror) for all optical elements
    # :param dens: a list with densities of o.e. materials. "?" is accepted for looking in the database
    # :param thick: a list with the thickness in mm for all o.e.'s. Only applicable for filters
    # :param angle: a list with the grazing angles in mrad for all o.e.'s. Only applicable for mirrors
    # :param roughness:a list with the roughness RMS in A for all o.e.'s. Only applicable for mirrors
    # :param output_file: name of the output file (default=None, no output file)
    # :return: a dictionary with the results
    # """
    optical_system_dictionary = xpower_calc(energies=e,source=spectral_power,
                    substance=["C",coating,coating],flags=[0,1,1],dens=[3.53,2.33,2.33],
                    thick=[thickness_diamond_mm,1,1],
                    angle=[0,incident_angle_mrad,incident_angle_mrad],roughness=[0,0,0],
                    output_file=None)

    for key in optical_system_dictionary.keys():
        print(key)

    print(optical_system_dictionary["info"])

    for i,ilabel in enumerate(optical_system_dictionary["labels"]):
        print(i,ilabel)

    # 0 Photon Energy [eV]
    # 1 Source
    # 2 [oe 1] Total CS cm2/g
    # 3 [oe 1] Mu cm^-1
    # 4 [oe 1] Transmitivity
    # 5 [oe 1] Absorption
    # 6 Intensity after oe #1
    # 7 [oe 2] 1-Re[n]=delta
    # 8 [oe 2] Im[n]=beta
    # 9 [oe 2] delta/beta
    # 10 [oe 2] Reflectivity-s
    # 11 [oe 2] Transmitivity
    # 12 Intensity after oe #2
    # 13 [oe 3] 1-Re[n]=delta
    # 14 [oe 3] Im[n]=beta
    # 15 [oe 3] delta/beta
    # 16 [oe 3] Reflectivity-s
    # 17 [oe 3] Transmitivity
    # 18 Intensity after oe #3

    print(optical_system_dictionary["data"].shape)

    # I would be interested in:
    #
    # - Total Power [W] emitted in the slit aperture: power_in_spectrum
    #
    # - Absorbed Power [W] by Diamond Window:  integral of col6-col1
    #
    # - Absorbed Power [W] for 1rst and 2nd mirrors: :  integral of col112-col6 and integral of col18-col12
    #
    # - Fitted parameters from the power density distribution calculated in a 5*5 mm slit aperture:
    #
    #     - Maximum value [W/mm2]
    #
    #     - Gaussian Fit parameters for both axis: FWHM [mm]


    I0 = numpy.trapz( optical_system_dictionary["data"][1,:], x=e, axis=-1)

    I1 = numpy.trapz( optical_system_dictionary["data"][6,:], x=e, axis=-1)

    I2 = numpy.trapz( optical_system_dictionary["data"][12,:], x=e, axis=-1)

    I3 = numpy.trapz( optical_system_dictionary["data"][18,:], x=e, axis=-1)


    print("Source power: ",I0)
    print(" after diamond: ",I1)
    print(" after M1: ",I2)
    print(" after M2: ",I3)

    out_dictionary["diamond_absorbed"] = I0-I1
    out_dictionary["m1_absorbed"] = I1-I2
    out_dictionary["m2_absorbed"] = I2-I3

    #
    # power density
    #

    h, v, p, code = xoppy_calc_undulator_power_density(ELECTRONENERGY=6.0,ELECTRONENERGYSPREAD=0.001,ELECTRONCURRENT=0.2,\
                                   ELECTRONBEAMSIZEH=SIGMAX,ELECTRONBEAMSIZEV=SIGMAY,\
                                   ELECTRONBEAMDIVERGENCEH=SIGMAXP,ELECTRONBEAMDIVERGENCEV=SIGMAYP,\
                                   PERIODID=undulator_period,NPERIODS=N,KV=K,DISTANCE=distance,GAPH=5e-3,GAPV=5e-3,\
                                   HSLITPOINTS=101,VSLITPOINTS=101,METHOD=2,USEEMITTANCES=1)

    if do_plot:
        plot_image(p,h,v,title="power density E = %d keV"%photon_energy)


    #
    # fit power density
    #

    print("============= Fitting power density to a 2D Gaussian. ==============\n")
    print("Please use these results with care: check if the original data looks like a Gaussian.")
    fit_parameters = fit_gaussian2d(p,h,v)
    print(info_params(fit_parameters))
    H,V = numpy.meshgrid(h,v)
    data_fitted = twoD_Gaussian( (H,V), *fit_parameters)
    power_in_spectrum = p.sum()*(h[1]-h[0])*(v[1]-v[0])
    print("  Total power in the calculated data [W]: ",power_in_spectrum)
    power_in_spectrum_fit = data_fitted.sum()*(h[1]-h[0])*(v[1]-v[0])
    print("  Total power in the fitted data [W]: ",power_in_spectrum_fit)
    # plot_image(data_fitted.reshape((h.size,v.size)),h, v,title="FIT")
    print("====================================================\n")

    if do_plot:
        data_fitted.shape = (h.size,v.size)
        plot_image(data_fitted,h,v,title="FITTED power density E = %d keV"%photon_energy)

    out_dictionary["fit_parameters"] = fit_parameters
    out_dictionary["fit_percent_difference"] = 100 * (power_in_spectrum_fit - power_in_spectrum) / power_in_spectrum

    return out_dictionary




if __name__ == "__main__":


    Energy_keV                    = [ 5   ,    7    ,     10   ,     12   ,     15   ,     20  ,      25   ,     30   ,     40  ]
    #Undulator                    = [ U32 ,    U27  ,     U27  ,     U27  ,     U27  ,     U27 ,      U32  ,     U32  ,     U32  ]
    lambda0_cm                    = [ 3.2 ,    2.7  ,     2.7  ,     2.7  ,     2.7  ,     2.7 ,      3.2  ,     3.2  ,     3.2  ]
    N                             = [ 72  ,    85   ,     85   ,     85   ,     85   ,     85  ,      72   ,     72   ,     72  ]
    K                             = [ 1.53,    1.292,     0.756,     1.846,     1.535,     1.12,      1.044,     1.274,     1.102  ]
    Diamond_window_thickness_mm   = [ 0.3 ,    0.3  ,     0.3  ,     0.3  ,     0.3  ,     0.3 ,      0.3  ,     0.3  ,     0.3  ]
    Distance_from_source_m        = [ 27  ,    27   ,     27   ,     27   ,     27   ,     27  ,      27   ,     27   ,     27  ]
    V_mm                          = [ 0.6 ,    0.6  ,     0.6  ,     0.6  ,     0.6  ,     0.6 ,      0.6  ,     0.6  ,     0.6  ]
    H_mm                          = [ 1.2 ,    1.2  ,     1    ,     0.9  ,     1.2  ,     1.2 ,      1.2  ,     0.9  ,     0.9  ]
    Coating                       = [ "Si"  ,  "Si"   ,   "Si"   ,   "Si"   ,   "Rh"   ,   "Rh"  ,    "Pt"   ,   "Pt"   ,   "Pt"  ]
    Incident_angle_mrad           = [ 4   ,    3    ,     2.5  ,     2    ,     4    ,     3   ,      3    ,     2    ,     2  ]


    #
    # calculation loop
    #
    out_dictionaries = []
    for i,photon_energy in enumerate(Energy_keV):
        out_dictionary = calculate_line(photon_energy,1e-2*lambda0_cm[i],N[i],K[i],Diamond_window_thickness_mm[i],
                                        Distance_from_source_m[i],1e-3*H_mm[i],1e-3*V_mm[i],Coating[i],Incident_angle_mrad[i],
                                        do_plot=False)
        out_dictionaries.append(out_dictionary)


    #
    # prepare text output
    #
    text_output = ""

    titles = ["energy_kev","power_in_spectrum","diamond_absorbed","m1_absorbed","m2_absorbed"]
    text_output += (" %20s %20s %20s %20s %20s \n")%(tuple(titles))

    for i in range(len(out_dictionaries)):
        text_output += ("%20d %20.3f %20.3f %20.3f %20.3f \n")%( Energy_keV[i],
                                                         out_dictionaries[i]["power_in_spectrum"],
                                                         out_dictionaries[i]["diamond_absorbed"],
                                                         out_dictionaries[i]["m1_absorbed"],
                                                         out_dictionaries[i]["m2_absorbed"])



    text_fit = ""
    titles_fit = ["energy_kev","Height A: ","center x0:","center y0","sigmax","sigmay","angle","offset","fit difference"]
    text_fit += ("%20s %20s %20s %20s %20s %20s %20s %20s %20s\n")%(tuple(titles_fit))
    for i in range(len(out_dictionaries)):
        text_fit += ("%20d %20.3f %20.3f %20.3f %20.3f %20.3f %20.3f %20.3f %20.3f \n")%(
                                                        Energy_keV[i],
                                                         out_dictionaries[i]["fit_parameters"][0],
                                                         out_dictionaries[i]["fit_parameters"][1],
                                                         out_dictionaries[i]["fit_parameters"][2],
                                                         out_dictionaries[i]["fit_parameters"][3],
                                                         out_dictionaries[i]["fit_parameters"][4],
                                                         out_dictionaries[i]["fit_parameters"][5],
                                                         out_dictionaries[i]["fit_parameters"][6],
                                                         out_dictionaries[i]["fit_percent_difference"])



    print(text_output)
    print(text_fit)

    #
    # dump to file
    #
    f = open("script1.txt",'w')
    f.write(text_output)
    f.write("\n\n\n")
    f.write(text_fit)
    f.close()

    print("File written to disk: script1.txt")
