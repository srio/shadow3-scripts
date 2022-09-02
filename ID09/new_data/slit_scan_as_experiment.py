#
# Python script to run shadow3 using libpyvinyl.
#
import numpy
from libpyvinyl.Parameters.Collections import CalculatorParameters
from libpyvinyl.Parameters.Parameter import Parameter
from shadow3libpyvinyl.Shadow3Calculator import Shadow3Calculator
from shadow3libpyvinyl.Shadow3Data import Shadow3BeamFormat, Shadow3OpenPMDFormat
from shadow3libpyvinyl.Shadow3Data import Shadow3Data



def add_source(parameters):
    #
    # Define variables. See https://raw.githubusercontent.com/oasys-kit/shadow3/master/docs/source.nml
    #

    p = Parameter('oe0.FDISTR',''); p.value = 3 ; parameters.add(p)
    p = Parameter('oe0.F_PHOT',''); p.value = 0 ; parameters.add(p)
    p = Parameter('oe0.HDIV1',''); p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe0.HDIV2',''); p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe0.IDO_VX',''); p.value = 0 ; parameters.add(p)
    p = Parameter('oe0.IDO_VZ',''); p.value = 0 ; parameters.add(p)
    p = Parameter('oe0.IDO_X_S',''); p.value = 0 ; parameters.add(p)
    p = Parameter('oe0.IDO_Y_S',''); p.value = 0 ; parameters.add(p)
    p = Parameter('oe0.IDO_Z_S',''); p.value = 0 ; parameters.add(p)
    p = Parameter('oe0.ISTAR1',''); p.value = 5676561 ; parameters.add(p)
    p = Parameter('oe0.NPOINT',''); p.value = 500000 ; parameters.add(p)
    p = Parameter('oe0.PH1',''); p.value = 19905.0 ; parameters.add(p)
    p = Parameter('oe0.SIGDIX',''); p.value = 1.75e-05 ; parameters.add(p)
    p = Parameter('oe0.SIGDIZ',''); p.value = 1.75e-05 ; parameters.add(p)
    p = Parameter('oe0.SIGMAX',''); p.value = 2.6e-05 ; parameters.add(p)
    p = Parameter('oe0.SIGMAZ',''); p.value = 6e-06 ; parameters.add(p)
    p = Parameter('oe0.VDIV1',''); p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe0.VDIV2',''); p.value = 0.0 ; parameters.add(p)

    return parameters
    

def add_beamline(parameters, use_slope_errors=True):

    
    #
    # Define variables. See https://raw.githubusercontent.com/oasys-kit/shadow3/master/docs/oe.nml
    #


    p = Parameter('oe1.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
    p = Parameter('oe1.FWRITE','') ; p.value = 3 ; parameters.add(p)
    p = Parameter('oe1.F_REFRAC','') ; p.value = 2 ; parameters.add(p)
    p = Parameter('oe1.F_SCREEN','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe1.I_SLIT','') ; p.value = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0]) ; parameters.add(p)
    p = Parameter('oe1.N_SCREEN','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe1.RX_SLIT','') ; p.value = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ; parameters.add(p)
    p = Parameter('oe1.RZ_SLIT','') ; p.value = numpy.array([0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ; parameters.add(p)
    p = Parameter('oe1.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe1.T_INCIDENCE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe1.T_REFLECTION','') ; p.value = 180.0 ; parameters.add(p)
    p = Parameter('oe1.T_SOURCE','') ; p.value = 27.066 ; parameters.add(p)

    p = Parameter('oe2.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
    p = Parameter('oe2.FILE_REFL','') ; p.value = '/nobackup/gurb1/srio/Oasys/Pd_reflec.dat' ; parameters.add(p)
    p = Parameter('oe2.FMIRR','') ; p.value = 3 ; parameters.add(p)
    p = Parameter('oe2.FWRITE','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe2.F_EXT','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe2.F_REFLEC','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe2.R_MAJ','') ; p.value = 20000.0 ; parameters.add(p)
    p = Parameter('oe2.R_MIN','') ; p.value = 0.046065 ; parameters.add(p)
    p = Parameter('oe2.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe2.T_INCIDENCE','') ; p.value = 89.8567605512 ; parameters.add(p)
    p = Parameter('oe2.T_REFLECTION','') ; p.value = 89.8567605512 ; parameters.add(p)
    p = Parameter('oe2.T_SOURCE','') ; p.value = 17.474 ; parameters.add(p)
    if use_slope_errors:
        p = Parameter('oe2.FHIT_C','') ; p.value = 1 ; parameters.add(p)
        p = Parameter('oe2.FILE_RIP','') ; p.value = '/nobackup/gurb1/srio/Oasys/mirror_shadow.dat' ; parameters.add(p)
        p = Parameter('oe2.F_G_S','') ; p.value = 2 ; parameters.add(p)
        p = Parameter('oe2.F_RIPPLE','') ; p.value = 1 ; parameters.add(p)
        p = Parameter('oe2.RLEN1','') ; p.value = 0.3 ; parameters.add(p)
        p = Parameter('oe2.RLEN2','') ; p.value = 0.3 ; parameters.add(p)
        p = Parameter('oe2.RWIDX1','') ; p.value = 0.065 ; parameters.add(p)
        p = Parameter('oe2.RWIDX2','') ; p.value = 0.065 ; parameters.add(p)

    p = Parameter('oe3.ALPHA','') ; p.value = 180.0 ; parameters.add(p)
    p = Parameter('oe3.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
    p = Parameter('oe3.FWRITE','') ; p.value = 3 ; parameters.add(p)
    p = Parameter('oe3.F_REFRAC','') ; p.value = 2 ; parameters.add(p)
    p = Parameter('oe3.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe3.T_INCIDENCE','') ; p.value = 180.0 ; parameters.add(p)
    p = Parameter('oe3.T_REFLECTION','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe3.T_SOURCE','') ; p.value = 0.0 ; parameters.add(p)

    p = Parameter('oe4.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
    p = Parameter('oe4.FWRITE','') ; p.value = 3 ; parameters.add(p)
    p = Parameter('oe4.F_REFRAC','') ; p.value = 2 ; parameters.add(p)
    p = Parameter('oe4.F_SCREEN','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe4.N_SCREEN','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe4.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe4.T_INCIDENCE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe4.T_REFLECTION','') ; p.value = 180.0 ; parameters.add(p)
    p = Parameter('oe4.T_SOURCE','') ; p.value = 9.819 ; parameters.add(p)

    # p = Parameter('oe5.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
    # p = Parameter('oe5.FWRITE','') ; p.value = 3 ; parameters.add(p)
    # p = Parameter('oe5.F_REFRAC','') ; p.value = 2 ; parameters.add(p)
    # p = Parameter('oe5.F_SCREEN','') ; p.value = 1 ; parameters.add(p)
    # p = Parameter('oe5.N_SCREEN','') ; p.value = 1 ; parameters.add(p)
    # p = Parameter('oe5.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
    # p = Parameter('oe5.T_INCIDENCE','') ; p.value = 0.0 ; parameters.add(p)
    # p = Parameter('oe5.T_REFLECTION','') ; p.value = 180.0 ; parameters.add(p)
    # p = Parameter('oe5.T_SOURCE','') ; p.value = 0.87 ; parameters.add(p)

    return parameters

def add_xeye(parameters):

    p = Parameter('oe5.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
    p = Parameter('oe5.FWRITE','') ; p.value = 3 ; parameters.add(p)
    p = Parameter('oe5.F_REFRAC','') ; p.value = 2 ; parameters.add(p)
    p = Parameter('oe5.F_SCREEN','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe5.N_SCREEN','') ; p.value = 1 ; parameters.add(p)
    p = Parameter('oe5.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe5.T_INCIDENCE','') ; p.value = 0.0 ; parameters.add(p)
    p = Parameter('oe5.T_REFLECTION','') ; p.value = 180.0 ; parameters.add(p)
    p = Parameter('oe5.T_SOURCE','') ; p.value = 0.87 ; parameters.add(p)

    return parameters

def run(two_steps=1, use_slope_errors=True):

    if two_steps:
        parameters = CalculatorParameters()
        parameters = add_source(parameters)
        parameters = add_beamline(parameters, use_slope_errors=use_slope_errors)

        calculator = Shadow3Calculator("", None, parameters=parameters)
        calculator.backengine()
        # calculator.data.write("tmp.h5", Shadow3OpenPMDFormat)  # openPMD data format
        output1 = calculator.output[calculator.output_keys[0]]


        # step two

        input_data =  output1.duplicate() # data.duplicate()

        parameters2 = CalculatorParameters()
        parameters2 = add_xeye(parameters2)
        print(parameters2)

        calculator2 = Shadow3Calculator("", input_data, parameters=parameters2)
        calculator2.backengine()

        output2 = calculator2.output[calculator2.output_keys[0]]

        return output1, output2

    else:
        parameters = CalculatorParameters()
        parameters = add_source(parameters)
        parameters = add_beamline(parameters, use_slope_errors=use_slope_errors)
        parameters = add_xeye(parameters)

        calculator = Shadow3Calculator("", None, parameters=parameters)
        calculator.backengine()

        output2 = calculator.output[calculator.output_keys[0]]

        return None, output2

#
# main
#
if __name__ == "__main__":
    import Shadow
    import os

    os.system("rm -f tmp.h5")
    from srxraylib.plot.gol import plot

    out = run(two_steps=1,use_slope_errors=True)

    if 0:
        if out[0] is not None:
            beam = Shadow.Beam(N=out[0].get_data()["nrays"])
            beam.rays = out[0].get_data()["rays"]
            Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space SLITTTT", )


        beam = Shadow.Beam(N=out[1].get_data()["nrays"])
        beam.rays = out[1].get_data()["rays"]
        Shadow.ShadowTools.plotxy(beam, 1, 3, nbins=101, nolost=1, title="Real space  XEYEEEEEE")


    #
    # analysis
    #

    do_plot = 0


    slit_aperture = 10e-6

    slit_positions = numpy.linspace(-600e-6, 600e-6, 300)
    xeye_positions = numpy.zeros_like(slit_positions)


    for i in range(len(slit_positions)):
        slit_position = slit_positions[i]

        # at slit
        # input_data = Shadow3Data("")
        # input_data.set_file("tmp.h5", Shadow3OpenPMDFormat)  # openPMD data format

        beamS = Shadow.Beam(N=out[0].get_data()["nrays"])
        beamS.rays = out[0].get_data()["rays"].copy()

        # at xeye
        beamX = Shadow.Beam(N=out[1].get_data()["nrays"])
        beamX.rays = out[1].get_data()["rays"].copy()


        z = beamS.getshonecol(3)
        bad = numpy.argwhere( numpy.abs(z-slit_position) > (slit_aperture / 2) )
        beamS.rays[bad,9] = -1001

        if do_plot: Shadow.ShadowTools.histo1(beamS, 3, nbins=1000, xrange=[-600e-6,600e-6], nolost=1)

        beamX.rays[:,9] = beamS.rays[:,9]
        tkt = beamX.histo1(3, nbins=1000, xrange=[-600e-6,600e-6], nolost=1)
        try:
            average = numpy.average(tkt['bin_center'], weights=tkt['histogram'])
        except:
            average = 0.0
        print("average position on Xeye: ", average)
        xeye_positions[i] = average

        if do_plot: plot(tkt['bin_center'],tkt['histogram'])

    plot(slit_positions, xeye_positions)

    fname = "slit_scan_as_experiment.dat"
    f = open(fname,"w")
    for i in range(slit_positions.size):
        f.write("%g  %g\n" % (slit_positions[i], xeye_positions[i]))
    f.close()
    print("file written to disk: %s" % fname)

