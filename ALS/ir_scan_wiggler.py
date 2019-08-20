#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy


def run_shadow(grazing=22.3):
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

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 1.9
    oe0.CONV_FACT = 1.0
    oe0.EPSI_X = 1.989e-09
    oe0.EPSI_Z = 3.007e-11
    oe0.FDISTR = 0
    oe0.FILE_TRAJ = b'/Users/srio/Oasys/xshwig.sha'
    oe0.FSOUR = 0
    oe0.FSOURCE_DEPTH = 0
    oe0.F_COLOR = 0
    oe0.F_PHOT = 0
    oe0.F_WIGGLER = 1
    oe0.HDIV1 = 1.0
    oe0.HDIV2 = 1.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 50000
    oe0.N_COLOR = 0
    oe0.PH1 = 1000.0
    oe0.PH2 = 1001.0
    oe0.POL_DEG = 0.0
    oe0.SIGMAX = 3.9e-05
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 3.1e-05
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.ALPHA = 90.0
    oe1.DUMMY = 100.0
    oe1.FCYL = 1
    oe1.FMIRR = 1
    oe1.FWRITE = 1
    oe1.F_DEFAULT = 0
    oe1.SIMAG = 4.29
    oe1.SSOUR = 1.58
    oe1.THETA = 90 - grazing
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 90 - grazing
    oe1.T_REFLECTION = 90 - grazing
    oe1.T_SOURCE = 1.58

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FCYL = 1
    oe2.FMIRR = 1
    oe2.FWRITE = 1
    oe2.F_DEFAULT = 0
    oe2.SIMAG = 3.11
    oe2.SSOUR = 2.76
    oe2.THETA = 45.0
    oe2.T_IMAGE = 3.11
    oe2.T_INCIDENCE = 45.0
    oe2.T_REFLECTION = 45.0
    oe2.T_SOURCE = 1.18

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

    # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")

    return beam, oe1

if __name__ == "__main__":

    from srxraylib.plot.gol import plot
    from srxraylib.util.h5_simple_writer import H5SimpleWriter

    Grazing = numpy.linspace(5,25,80)
    # Grazing = numpy.linspace(15, 88, 50)
    Radius = numpy.zeros_like(Grazing)
    Fwhm = numpy.zeros_like(Grazing)
    Std = numpy.zeros_like(Grazing)

    h = H5SimpleWriter.initialize_file("ir_scan_wiggler.h5")

    for i,grazing in enumerate(Grazing):

        beam,oe1 = run_shadow(grazing=grazing)
        # Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space")

        tkt = beam.histo1(1,nbins=201,nolost=1) # xrange=[-4000e-6,4000e-6],
        print(grazing,oe1.RMIRR,1e6*tkt["fwhm"])
        Radius[i] = oe1.RMIRR
        Fwhm[i] = 1e6*tkt["fwhm_subpixel"]
        Std[i] = 1e6 * beam.get_standard_deviation(1,nolost=1)

        h.create_entry("iteration %f" % grazing, nx_default="histogram")
        h.add_dataset(1e6*tkt["bin_center"],tkt["histogram"], dataset_name="histogram",entry_name="iteration %f" % grazing,
                      title_x="X / um",title_y="intensity / a.u.")



    plot(Radius,Fwhm,xtitle="R/m",ytitle="FWHM/um",show=False)
    plot(Grazing, Fwhm,Grazing, Std,xtitle="Graz/deg",ytitle="FWHM/um",legend=["fwhm","std"])

    h.create_entry("scan results", nx_default="Fwhm")
    h.add_dataset(Grazing, Fwhm, dataset_name="Fwhm", entry_name="scan results",
                  title_x="grazing angle / deg", title_y="fwhm / um")
    h.add_dataset(Grazing, Std, dataset_name="Std", entry_name="scan results",
                  title_x="grazing angle / deg", title_y="fwhm / um")
    h.add_dataset(Radius, Fwhm, dataset_name="Radius", entry_name="scan results",
                  title_x="radius / m", title_y="fwhm / um")

    for key in tkt.keys():
        print(key)
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
