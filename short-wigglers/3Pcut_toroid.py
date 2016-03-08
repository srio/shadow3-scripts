#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy

# write (1) or not (0) SHADOW files start.xx end.xx star.xx
iwrite = 1

#
# initialize shadow3 source (oe0) and beam
#
beam = Shadow.Beam()
oe0 = Shadow.Source()
oe1 = Shadow.OE()

#
# Define variables. See meaning of variables in: 
#  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
#  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
#

oe0.BENER = 6.0
oe0.CONV_FACT = 100.0
oe0.EPSI_DX = 89.4
oe0.EPSI_DZ = -104.8
oe0.EPSI_X = 2.16e-08
oe0.EPSI_Z = 5e-10
oe0.FDISTR = 0
oe0.FILE_TRAJ = b'/Users/srio/Oasys/OasysRun/xshwig.sha'
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
oe0.NPOINT = 20000
oe0.NTOTALPOINT = 0
oe0.N_COLOR = 0
oe0.PH1 = 20000.0
oe0.PH2 = 20040.0
oe0.POL_DEG = 0.0
oe0.SIGMAX = 0.0008757
oe0.SIGMAY = 0.0
oe0.SIGMAZ = 0.0001647
oe0.VDIV1 = 1.0
oe0.VDIV2 = 1.0
oe0.WXSOU = 0.0
oe0.WYSOU = 0.0
oe0.WZSOU = 0.0

oe1.DUMMY = 1.0
oe1.FHIT_C = 1
oe1.FMIRR = 3
oe1.FWRITE = 1
oe1.RLEN1 = 50.0
oe1.RLEN2 = 50.0
oe1.RWIDX1 = 10.0
oe1.RWIDX2 = 10.0
oe1.T_IMAGE = 3000.0
oe1.T_INCIDENCE = 89.828
oe1.T_REFLECTION = 89.828
oe1.T_SOURCE = 3000.0



#Run SHADOW to create the source

if iwrite:
    oe0.write("start.00")

beam.genSource(oe0)

if iwrite:
    oe0.write("end.00")
    beam.write("begin.dat")


#
#run optical element 1
#
print("    Running optical element: %d"%(1))
if iwrite:
    oe1.write("start.01")
beam.traceOE(oe1,1)
if iwrite:
    oe1.write("end.01")
    beam.write("star.01")


Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
# Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
# Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    
