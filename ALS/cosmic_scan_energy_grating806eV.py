import Shadow
import numpy
from source_tools import get_sigmas_ALSU

def run_shadow(energy=806.0, delta_energy=0.0, alpha = 88.603245, beta = 87.796714, sigmas=[1e-5,1e-5,1e-5]):
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

if __name__ == "__main__":

    from srxraylib.plot.gol import plot

    from grating_tools import trajectories
    from respower import respower, respower_plot


    energy1 = 200.0 # 5000.0
    undulator_length = 2.1

    Energy = numpy.linspace(250.0, 2500.0, 20)
    Resolution = Energy * 0.0
    Alpha = Energy * 0.0
    Beta = Energy * 0.0
    ImageSize = Energy * 0.0
    SourceSize = Energy * 0.0
    M = Energy * 0.0
    C = Energy * 0.0


    for i,energy1 in enumerate(Energy):
        # trajectories(energies, r, rp, k0, m, b2, verbose=False)

        sigmas = get_sigmas_ALSU(energy1,undulator_length)
        alpha1,beta1 = trajectories(energy1, 25.201, 7.573, 287444, 1, b2=0.24736325719129112, verbose=False)
        alpha1 *= 180.0 / numpy.pi
        beta1 *= -180.0 / numpy.pi


        #
        # print("alpha: ",alpha,alpha1)
        # print("beta: ",beta, beta1)


        beam = run_shadow(energy1,0.4,alpha1,beta1,sigmas)
        dict = respower(beam, 11, 1, nolost=True)

        print("\n\nset energy: %f eV"%energy1)
        print("alpha: ",alpha1)
        print("beta: ",beta1)
        print("Resolving power: %f"%dict["resolvingPower"])


        # store results
        Energy[i] = energy1
        Resolution[i] = dict["resolvingPower"]
        Alpha[i] = alpha1
        Beta[i] = beta1

        # run now the monochromatic case to get the monochromatic size
        beam1 = run_shadow(energy1, 0.0,  alpha1, beta1, sigmas)
        # tkt = dict["histo_dict"]
        tkt = beam1.histo1(1,nolost=True, ref=23, nbins=100)
        # plot(tkt["bin_path"] * 1e6, tkt["histogram_path"], title="E=%s  S=%s"%(energy1,tkt["fwhm"] * 1e6))
        ImageSize[i] = tkt["fwhm"]*1e6

        SourceSize[i] = 2.35 * sigmas[1] * 1e6

        C[i] = numpy.cos(beta1*numpy.pi/180) / numpy.cos(alpha1*numpy.pi/180)
        M[i] = 7.573 / 25.201 / C[i]




    # plot(Energy,Resolution)

    # dump to file
    dumpfile = "cosmic_scan_energy_grating806eV.dat"
    f = open(dumpfile,'w')
    f.write("\n#S 1\n")
    f.write("#N 8\n")
    f.write("#L photon energy [eV]  alpha [deg]  beta [deg]  size source V [um]  size [um]  source*M  C  resolving power\n")
    for i in range(Energy.size):
        f.write("%f   %f   %f   %f   %f   %f   %f   %f\n"%(Energy[i],Alpha[i],Beta[i],SourceSize[i],
                                                           ImageSize[i],SourceSize[i]*M[i],C[i],Resolution[i]))
    f.close()
    print("File written to disk: %s"%dumpfile)