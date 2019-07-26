import Shadow
import numpy
from source_tools import get_sigmas_ALSU

##############################################################################################################
##############################################################################################################
##############################################################################################################

def run_shadow_E379(energy=806.0, delta_energy=0.0, alpha = 88.603245, beta = 87.796714, sigmas=[1e-5,1e-5,1e-5]):
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 50000
    oe0.PH1 = 379.1
    oe0.PH2 = 379.1
    oe0.SIGDIX = 2.8130065349342457e-05
    oe0.SIGDIZ = 2.7944598343119502e-05
    oe0.SIGMAX = 2.1574735780763598e-05
    oe0.SIGMAZ = 2.3133292545804222e-05
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.ALPHA = 270.0
    oe1.DUMMY = 100.0
    oe1.FHIT_C = 1
    oe1.FWRITE = 1
    oe1.RLEN1 = 0.05
    oe1.RLEN2 = 0.05
    oe1.RWIDX1 = 0.005
    oe1.RWIDX2 = 0.005
    oe1.T_IMAGE = -15.098
    oe1.T_INCIDENCE = 88.75
    oe1.T_REFLECTION = 88.75
    oe1.T_SOURCE = 15.098

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FWRITE = 1
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 87.999955
    oe2.T_REFLECTION = 87.999955
    oe2.T_SOURCE = 25.201

    oe3.ALPHA = 180.0
    oe3.DUMMY = 100.0
    oe3.FWRITE = 1
    oe3.F_GRATING = 1
    oe3.F_RULING = 5
    oe3.F_RUL_ABS = 1
    oe3.RULING = 178960.0
    oe3.RUL_A1 = 84159.42650645
    oe3.RUL_A2 = 14457.8147684
    oe3.RUL_A3 = 2658.03276251
    oe3.T_IMAGE = 7.573
    oe3.T_INCIDENCE = 88.480393
    oe3.T_REFLECTION = 87.519517
    oe3.T_SOURCE = 0.0

    oe4.ALPHA = 90.0
    oe4.DUMMY = 100.0
    oe4.FCYL = 1
    oe4.FMIRR = 2
    oe4.FWRITE = 1
    oe4.F_DEFAULT = 0
    oe4.SIMAG = 5.573
    oe4.SSOUR = 27.201
    oe4.THETA = 88.75
    oe4.T_IMAGE = 5.573
    oe4.T_INCIDENCE = 88.75
    oe4.T_REFLECTION = 88.75
    oe4.T_SOURCE = -5.573

    oe0.NPOINT = 500000
    oe0.PH1 = energy - 0.5*delta_energy
    oe0.PH2 = energy + 0.5*delta_energy
    oe0.SIGDIX = sigmas[2] # 1.9733083328114146e-05
    oe0.SIGDIZ = sigmas[3] # 1.946778306932498e-05
    oe0.SIGMAX = sigmas[0] # 1.7218556115924652e-05
    oe0.SIGMAZ = sigmas[1] # 1.913527305050143e-05
    oe2.T_INCIDENCE = 0.5*(alpha+beta) # 88.1999795
    oe2.T_REFLECTION = 0.5*(alpha+beta) #88.1999795
    oe3.T_INCIDENCE = alpha # 88.603245
    oe3.T_REFLECTION = beta # 87.796714
    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    #
    # run optical element 4
    #
    print("    Running optical element: %d" % (4))
    if iwrite:
        oe4.write("start.04")

    beam.traceOE(oe4, 4)

    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return beam

##############################################################################################################
##############################################################################################################
##############################################################################################################

def run_shadow_E1714(energy=806.0, delta_energy=0.0, alpha=88.603245, beta=87.796714, sigmas=[1e-5, 1e-5, 1e-5]):

    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 50000
    oe0.PH1 = 1714.4
    oe0.PH2 = 1714.4
    oe0.SIGDIX = 1.4152153523560814e-05
    oe0.SIGDIZ = 1.3779820367277461e-05
    oe0.SIGMAX = 1.4729649862003558e-05
    oe0.SIGMAZ = 1.6930522291329982e-05
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.ALPHA = 270.0
    oe1.DUMMY = 100.0
    oe1.FHIT_C = 1
    oe1.FWRITE = 1
    oe1.RLEN1 = 0.05
    oe1.RLEN2 = 0.05
    oe1.RWIDX1 = 0.005
    oe1.RWIDX2 = 0.005
    oe1.T_IMAGE = -15.098
    oe1.T_INCIDENCE = 88.75
    oe1.T_REFLECTION = 88.75
    oe1.T_SOURCE = 15.098

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FWRITE = 1
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 88.750009
    oe2.T_REFLECTION = 88.750009
    oe2.T_SOURCE = 25.201

    oe3.ALPHA = 180.0
    oe3.DUMMY = 100.0
    oe3.FWRITE = 1
    oe3.F_GRATING = 1
    oe3.F_RULING = 5
    oe3.F_RUL_ABS = 1
    oe3.RULING = 352450.0
    oe3.RUL_A1 = 153646.17916295
    oe3.RUL_A2 = 26816.18910246
    oe3.RUL_A3 = 4909.34710357
    oe3.T_IMAGE = 7.573
    oe3.T_INCIDENCE = 89.084741
    oe3.T_REFLECTION = 88.415277
    oe3.T_SOURCE = 0.0

    oe4.ALPHA = 90.0
    oe4.DUMMY = 100.0
    oe4.FCYL = 1
    oe4.FMIRR = 2
    oe4.FWRITE = 1
    oe4.F_DEFAULT = 0
    oe4.SIMAG = 5.573
    oe4.SSOUR = 27.201
    oe4.THETA = 88.75
    oe4.T_IMAGE = 5.573
    oe4.T_INCIDENCE = 88.75
    oe4.T_REFLECTION = 88.75
    oe4.T_SOURCE = -5.573

    oe0.NPOINT = 500000
    oe0.PH1 = energy - 0.5*delta_energy
    oe0.PH2 = energy + 0.5*delta_energy
    oe0.SIGDIX = sigmas[2] # 1.9733083328114146e-05
    oe0.SIGDIZ = sigmas[3] # 1.946778306932498e-05
    oe0.SIGMAX = sigmas[0] # 1.7218556115924652e-05
    oe0.SIGMAZ = sigmas[1] # 1.913527305050143e-05
    oe2.T_INCIDENCE = 0.5*(alpha+beta) # 88.1999795
    oe2.T_REFLECTION = 0.5*(alpha+beta) #88.1999795
    oe3.T_INCIDENCE = alpha # 88.603245
    oe3.T_REFLECTION = beta # 87.796714

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    #
    # run optical element 4
    #
    print("    Running optical element: %d" % (4))
    if iwrite:
        oe4.write("start.04")

    beam.traceOE(oe4, 4)

    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    return beam

def run_shadow_E806(energy=806.0, delta_energy=0.0, alpha = 88.603245, beta = 87.796714, sigmas=[1e-5,1e-5,1e-5]):
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #


    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 500000
    oe0.PH1 = energy - 0.5*delta_energy
    oe0.PH2 = energy + 0.5*delta_energy
    oe0.SIGDIX = sigmas[2] # 1.9733083328114146e-05
    oe0.SIGDIZ = sigmas[3] # 1.946778306932498e-05
    oe0.SIGMAX = sigmas[0] # 1.7218556115924652e-05
    oe0.SIGMAZ = sigmas[1] # 1.913527305050143e-05
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.ALPHA = 270.0
    oe1.DUMMY = 100.0
    oe1.FHIT_C = 1
    oe1.FWRITE = 1
    oe1.RLEN1 = 0.05
    oe1.RLEN2 = 0.05
    oe1.RWIDX1 = 0.005
    oe1.RWIDX2 = 0.005
    oe1.T_IMAGE = -15.098
    oe1.T_INCIDENCE = 88.75
    oe1.T_REFLECTION = 88.75
    oe1.T_SOURCE = 15.098

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FWRITE = 1
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 0.5*(alpha+beta) # 88.1999795
    oe2.T_REFLECTION = 0.5*(alpha+beta) #88.1999795
    oe2.T_SOURCE = 25.201

    oe3.ALPHA = 180.0
    oe3.DUMMY = 100.0
    oe3.FWRITE = 1
    oe3.F_GRATING = 1
    oe3.F_RULING = 5
    oe3.F_RUL_ABS = 1
    oe3.RULING = 287440.0
    oe3.RUL_A1 = 142204.18929413
    oe3.RUL_A2 = 24200.21342526
    oe3.RUL_A3 = 4464.69044186
    oe3.T_IMAGE = 7.573
    oe3.T_INCIDENCE = alpha # 88.603245
    oe3.T_REFLECTION = beta # 87.796714
    oe3.T_SOURCE = 0.0

    oe4.ALPHA = 90.0
    oe4.DUMMY = 100.0
    oe4.FCYL = 1
    oe4.FMIRR = 2
    oe4.FWRITE = 1
    oe4.F_DEFAULT = 0
    oe4.SIMAG = 5.573
    oe4.SSOUR = 27.201
    oe4.THETA = 88.75
    oe4.T_IMAGE = 5.573
    oe4.T_INCIDENCE = 88.75
    oe4.T_REFLECTION = 88.75
    oe4.T_SOURCE = -5.573

    #
    # overwrite varibles
    #
    oe0.NPOINT = 500000
    oe0.PH1 = energy - 0.5*delta_energy
    oe0.PH2 = energy + 0.5*delta_energy
    oe0.SIGDIX = sigmas[2] # 1.9733083328114146e-05
    oe0.SIGDIZ = sigmas[3] # 1.946778306932498e-05
    oe0.SIGMAX = sigmas[0] # 1.7218556115924652e-05
    oe0.SIGMAZ = sigmas[1] # 1.913527305050143e-05
    oe2.T_INCIDENCE = 0.5*(alpha+beta) # 88.1999795
    oe2.T_REFLECTION = 0.5*(alpha+beta) #88.1999795
    oe3.T_INCIDENCE = alpha # 88.603245
    oe3.T_REFLECTION = beta # 87.796714


    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    #
    # run optical element 2
    #
    print("    Running optical element: %d" % (2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2, 2)

    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")

    #
    # run optical element 3
    #
    print("    Running optical element: %d" % (3))
    if iwrite:
        oe3.write("start.03")

    beam.traceOE(oe3, 3)

    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")

    #
    # run optical element 4
    #
    print("    Running optical element: %d" % (4))
    if iwrite:
        oe4.write("start.04")

    beam.traceOE(oe4, 4)

    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")


    return beam

##############################################################################################################
##############################################################################################################
##############################################################################################################

if __name__ == "__main__":

    from srxraylib.plot.gol import plot

    from grating_tools import trajectories
    from respower import respower, respower_plot


    grating = 2

    undulator_length = 2.1

    Energy = numpy.linspace(250.0, 2500.0, 20)
    Resolution = Energy * 0.0
    Resolution5 = Energy * 0.0
    Alpha = Energy * 0.0
    Beta = Energy * 0.0
    ImageSize = Energy * 0.0
    SourceSize = Energy * 0.0
    Footprint = Energy * 0.0
    M = Energy * 0.0
    C = Energy * 0.0


    if grating == 1:
        run_shadow = run_shadow_E379
        dumpfile = "cosmic_scan_energy_grating379eV.dat"
        k0 = 178960.0
        b2 = 0.2351347410216088
    elif grating == 2:
        run_shadow = run_shadow_E806
        dumpfile = "cosmic_scan_energy_grating806eV.dat"
        k0 = 287444
        b2 = 0.24736325719129112
    elif grating == 3:
        run_shadow = run_shadow_E1714
        dumpfile = "cosmic_scan_energy_grating1714eV.dat"
        k0 = 352450.0
        b2 = 0.21796876033899015

    for i,energy1 in enumerate(Energy):
        # trajectories(energies, r, rp, k0, m, b2, verbose=False)

        sigmas = get_sigmas_ALSU(energy1,undulator_length)
        alpha1,beta1 = trajectories(energy1, 25.201, 7.573, k0, 1, b2=b2, verbose=False)
        alpha1 *= 180.0 / numpy.pi
        beta1 *= -180.0 / numpy.pi


        beam = run_shadow(energy1,0.4,alpha1,beta1,sigmas)

        print("\n\nset energy: %f eV"%energy1)
        print("alpha: ",alpha1)
        print("beta: ",beta1)

        dict = respower(beam, 11, 1, hlimit=0.1, nolost=True)
        dict5 = respower(beam, 11, 1, hlimit=0.5, nolost=True)


        # store results
        Energy[i] = energy1
        Resolution[i] = dict["resolvingPower"]
        Resolution5[i] = dict5["resolvingPower"]
        Alpha[i] = alpha1
        Beta[i] = beta1
        C[i] = numpy.cos(beta1*numpy.pi/180) / numpy.cos(alpha1*numpy.pi/180)
        M[i] = 7.573 / 25.201 / C[i]


        # run now the monochromatic case to get the monochromatic size
        # beam1 = run_shadow(energy1, 0.0,  alpha1, beta1, sigmas)
        # tkt = beam1.histo1(1,nolost=True, ref=23, nbins=100)
        tkt = dict["histo_dict"]
        # plot(tkt["bin_path"] * 1e6, tkt["histogram_path"], title="E=%s  S=%s"%(energy1,tkt["fwhm"] * 1e6))
        ImageSize[i] = tkt["fwhm"]*1e6
        SourceSize[i] = 2.35 * sigmas[1] * 1e6

        # get now the footprint on grating
        b = Shadow.Beam()
        b.load("mirr.03")
        tkt_m = b.histo1(2, nolost=True)
        Footprint[i] = tkt_m["fwhm"]


    #
    # dump to file
    #

    f = open(dumpfile,'w')
    f.write("\n#S 1\n")
    f.write("#N 10\n")
    f.write("#L photon energy [eV]  alpha [deg]  beta [deg]  size source V [um]  size [um]  source*M  footprint on grating [m]  C  resolving power 0.1  resolving power 0.5\n")
    for i in range(Energy.size):
        f.write("%f   %f   %f   %f   %f   %f   %f   %f   %f   %f\n"%(
            Energy[i],
            Alpha[i],
            Beta[i],
            SourceSize[i],
            ImageSize[i],
            SourceSize[i]*M[i],
            Footprint[i],
            C[i],
            Resolution[i],
            Resolution5[i]))
    f.close()
    print("File written to disk: %s"%dumpfile)