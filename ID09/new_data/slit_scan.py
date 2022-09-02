#
# Python script to run shadow3 using libpyvinyl.
#
import numpy
from libpyvinyl.Parameters.Collections import CalculatorParameters
from libpyvinyl.Parameters.Parameter import Parameter
from shadow3libpyvinyl.Shadow3Calculator import Shadow3Calculator
from shadow3libpyvinyl.Shadow3Data import Shadow3BeamFormat, Shadow3OpenPMDFormat



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
    

def add_beamline(parameters, use_slope_errors=False):

    
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



    if use_slope_errors:
        p = Parameter('oe2.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
        p = Parameter('oe2.FHIT_C','') ; p.value = 1 ; parameters.add(p)
        p = Parameter('oe2.FILE_REFL','') ; p.value = '/nobackup/gurb1/srio/Oasys/Pd_reflec.dat' ; parameters.add(p)
        p = Parameter('oe2.FILE_RIP','') ; p.value = '/nobackup/gurb1/srio/Oasys/mirror_shadow.dat' ; parameters.add(p)
        p = Parameter('oe2.FMIRR','') ; p.value = 3 ; parameters.add(p)
        p = Parameter('oe2.FWRITE','') ; p.value = 1 ; parameters.add(p)
        p = Parameter('oe2.F_EXT','') ; p.value = 1 ; parameters.add(p)
        p = Parameter('oe2.F_G_S','') ; p.value = 2 ; parameters.add(p)
        p = Parameter('oe2.F_REFLEC','') ; p.value = 1 ; parameters.add(p)
        p = Parameter('oe2.F_RIPPLE','') ; p.value = 1 ; parameters.add(p)
        p = Parameter('oe2.RLEN1','') ; p.value = 0.3 ; parameters.add(p)
        p = Parameter('oe2.RLEN2','') ; p.value = 0.3 ; parameters.add(p)
        p = Parameter('oe2.RWIDX1','') ; p.value = 0.065 ; parameters.add(p)
        p = Parameter('oe2.RWIDX2','') ; p.value = 0.065 ; parameters.add(p)
        p = Parameter('oe2.R_MAJ','') ; p.value = 20000.0 ; parameters.add(p)
        p = Parameter('oe2.R_MIN','') ; p.value = 0.046065 ; parameters.add(p)
        p = Parameter('oe2.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
        p = Parameter('oe2.T_INCIDENCE','') ; p.value = 89.8567605512 ; parameters.add(p)
        p = Parameter('oe2.T_REFLECTION','') ; p.value = 89.8567605512 ; parameters.add(p)
        p = Parameter('oe2.T_SOURCE','') ; p.value = 17.474 ; parameters.add(p)
    else:
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
    p = Parameter('oe4.T_SOURCE','') ; p.value = 10.689 ; parameters.add(p)

    # p = Parameter('oe5.CZ_SLIT','') ; p.value = numpy.array([0.0002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ; parameters.add(p)
    # p = Parameter('oe5.DUMMY','') ; p.value = 100.0 ; parameters.add(p)
    # p = Parameter('oe5.FWRITE','') ; p.value = 3 ; parameters.add(p)
    # p = Parameter('oe5.F_REFRAC','') ; p.value = 2 ; parameters.add(p)
    # p = Parameter('oe5.F_SCREEN','') ; p.value = 1 ; parameters.add(p)
    # p = Parameter('oe5.I_SLIT','') ; p.value = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0]) ; parameters.add(p)
    # p = Parameter('oe5.N_SCREEN','') ; p.value = 1 ; parameters.add(p)
    # p = Parameter('oe5.RX_SLIT','') ; p.value = numpy.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ; parameters.add(p)
    # p = Parameter('oe5.RZ_SLIT','') ; p.value = numpy.array([1e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) ; parameters.add(p)
    # p = Parameter('oe5.T_IMAGE','') ; p.value = 0.0 ; parameters.add(p)
    # p = Parameter('oe5.T_INCIDENCE','') ; p.value = 0.0 ; parameters.add(p)
    # p = Parameter('oe5.T_REFLECTION','') ; p.value = 180.0 ; parameters.add(p)
    # p = Parameter('oe5.T_SOURCE','') ; p.value = 0.0 ; parameters.add(p)



    return parameters

    
#
# main
#
import os
parameters = CalculatorParameters()
parameters = add_source(parameters)
parameters = add_beamline(parameters, use_slope_errors=False)

calculator  = Shadow3Calculator("", None, parameters=parameters)
calculator.backengine()
os.system("cp mirr.02 mirrNo.02")

parameters2 = CalculatorParameters()
parameters2 = add_source(parameters)
parameters2 = add_beamline(parameters, use_slope_errors=True)

calculator2  = Shadow3Calculator("", None, parameters=parameters2)
calculator2.backengine()

#
# output files 
#
# calculator.parameters.to_json("my_parameters.json")
# print(calculator.parameters)

# calculator.data.write("tmp.dat", Shadow3BeamFormat)    # raw data format
# calculator.data.write("tmp.h5", Shadow3OpenPMDFormat)  # openPMD data format

    
    
#
# analysis
#
import Shadow
from srxraylib.plot.gol import plot

do_plot = 0


slit_aperture = 10e-6

slit_positions = numpy.linspace(-600e-6, 600e-6, 300)
mirror_positions = numpy.zeros_like(slit_positions)
mirror_positions2 = numpy.zeros_like(slit_positions)


# calculator.data.write("star.04", Shadow3BeamFormat)

for i in range(len(slit_positions)):
    # slit_position = -300e-6
    slit_position = slit_positions[i]

    beam = Shadow.Beam(N=calculator.data.get_data()["nrays"])
    beam.rays = calculator.data.get_data()["rays"].copy()

    beam2 = Shadow.Beam(N=calculator.data.get_data()["nrays"])
    beam2.rays = calculator2.data.get_data()["rays"].copy()

    # if do_plot: Shadow.ShadowTools.histo1(beam, 3, nbins=1000, xrange=[-600e-6,600e-6], nolost=1)
    #reflag

    z = beam.getshonecol(3)
    bad = numpy.argwhere( numpy.abs(z-slit_position) > (slit_aperture / 2) )
    beam.rays[bad,9] = -1001

    z2 = beam2.getshonecol(3)
    bad = numpy.argwhere( numpy.abs(z2-slit_position) > (slit_aperture / 2) )
    beam2.rays[bad,9] = -1001

    if do_plot: Shadow.ShadowTools.histo1(beam, 3, nbins=1000, xrange=[-600e-6,600e-6], nolost=1)

    m = Shadow.Beam()
    m.load("mirrNo.02")
    m.rays[:,9] = beam.rays[:,9] # copy flag
    tkt = m.histo1(2, nbins=1000, xrange=[-0.4,0.4], nolost=1)
    try:
        average = numpy.average(tkt['bin_center'], weights=tkt['histogram'])
    except:
        average = 0.0
    print("average position on mirror (no slope): ", average)
    mirror_positions[i] = average

    m2 = Shadow.Beam()
    m2.load("mirr.02")
    m2.rays[:, 9] = beam2.rays[:, 9]  # copy flag
    tkt2 = m2.histo1(2, nbins=1000, xrange=[-0.4, 0.4], nolost=1)
    try:
        average = numpy.average(tkt2['bin_center'], weights=tkt2['histogram'])
    except:
        average = 0.0
    print("i=%d\naverage position on mirror: %f" % (i,average))
    mirror_positions2[i] = average

    if do_plot: plot(tkt['bin_center'],tkt['histogram'])

plot(slit_positions, mirror_positions,
     slit_positions, mirror_positions2)

fname = "slit_scan.dat"
f = open(fname,"w")
for i in range(slit_positions.size):
    f.write("%g  %g  %g\n" % (slit_positions[i], mirror_positions[i], mirror_positions2[i]))
f.close()
print("file written to disk: %s" % fname)