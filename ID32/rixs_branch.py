#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
import matplotlib.pyplot as plt

import copy

def trace_rixs_branch(fname,theta_mrad=19.1):
    
    
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
    oe5 = Shadow.OE()
    oe6 = Shadow.OE()
    oe7 = Shadow.OE()
    oe8 = Shadow.OE()
    oe9 = Shadow.OE()
    oe10 = Shadow.OE()
    
    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #
    
    # Source
    oe0.FDISTR = 3
    oe0.FSOURCE_DEPTH = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.F_POLAR = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 9553
    oe0.NPOINT = 500000
    oe0.PH1 = 929.95
    oe0.PH2 = 930.05
    oe0.SIGDIX = 1.66e-05
    oe0.SIGDIZ = 1.61e-05
    oe0.SIGMAX = 0.033
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 0.013
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0
    
    # Primary slits
    oe1.DUMMY = 0.1
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.RX_SLIT = numpy.array([3.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.RZ_SLIT = numpy.array([2.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.T_IMAGE = 1337.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 27275.0
    
    # Plane mirror
    oe2.ALPHA = 90.0
    oe2.DUMMY = 0.1
    oe2.FHIT_C = 1
    if fname:
        oe2.FILE_RIP = fname
        oe2.F_G_S = 2
        oe2.F_RIPPLE = 1
    oe2.RLEN1 = 130.0
    oe2.RLEN2 = 130.0
    oe2.RWIDX1 = 30.0
    oe2.RWIDX2 = 30.0
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 89.25
    oe2.T_REFLECTION = 89.25
    oe2.T_SOURCE = 1337.0
    
    # Toroidal mirror
    oe3.ALPHA = 180.0
    oe3.DUMMY = 0.1
    oe3.FHIT_C = 1
    oe3.FMIRR = 3
    oe3.F_EXT = 1
    oe3.RLEN1 = 130.0
    oe3.RLEN2 = 130.0
    oe3.RWIDX1 = 10.0
    oe3.RWIDX2 = 10.0
    oe3.R_MAJ = 4811761.2088
    oe3.R_MIN = 217.9
    oe3.T_IMAGE = 937.5
    oe3.T_INCIDENCE = 89.25
    oe3.T_REFLECTION = 89.25
    oe3.T_SOURCE = 611.0
    
    # Secondary slits
    oe4.DUMMY = 0.1
    oe4.FWRITE = 3
    oe4.F_REFRAC = 2
    oe4.F_SCREEN = 1
    oe4.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe4.N_SCREEN = 1
    oe4.RX_SLIT = numpy.array([3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe4.RZ_SLIT = numpy.array([3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe4.T_IMAGE = 4782.5
    oe4.T_INCIDENCE = 0.0
    oe4.T_REFLECTION = 180.0
    oe4.T_SOURCE = 937.5
    
    # Entrance slit
    oe5.DUMMY = 0.1
    oe5.FWRITE = 3
    oe5.F_REFRAC = 2
    oe5.F_SCREEN = 1
    oe5.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe5.N_SCREEN = 1
    oe5.RX_SLIT = numpy.array([0.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe5.RZ_SLIT = numpy.array([20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe5.T_IMAGE = 0.0
    oe5.T_INCIDENCE = 0.0
    oe5.T_REFLECTION = 180.0
    oe5.T_SOURCE = 4782.5
    
    # PGM plane mirror
    oe6.ALPHA = 270.0
    oe6.DUMMY = 0.1
    oe6.T_IMAGE = 116.3105
    oe6.T_INCIDENCE = 88.151426
    oe6.T_REFLECTION = 88.151426
    oe6.T_SOURCE = 22267.863
    
    # PGM VLS grating
    oe7.ALPHA = 180.0
    oe7.DUMMY = 0.1
    oe7.FHIT_C = 1
    oe7.F_GRATING = 1
    oe7.F_RULING = 5
    oe7.F_RUL_ABS = 1
    oe7.RLEN1 = 90.0
    oe7.RLEN2 = 90.0
    oe7.RULING = 800.0
    oe7.RUL_A1 = 0.05923758
    oe7.RUL_A2 = 1.63e-06
    oe7.RWIDX1 = 10.0
    oe7.RWIDX2 = 10.0
    oe7.T_IMAGE = 35000.0
    oe7.T_INCIDENCE = 89.098637
    oe7.T_REFLECTION = 87.204215
    oe7.T_SOURCE = 116.3105
    
    # Exit slit
    oe8.DUMMY = 0.1
    oe8.FWRITE = 3
    oe8.F_REFRAC = 2
    oe8.F_SCREEN = 1
    oe8.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe8.N_SCREEN = 1
    oe8.RX_SLIT = numpy.array([20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe8.RZ_SLIT = numpy.array([0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe8.T_IMAGE = 0.0
    oe8.T_INCIDENCE = 0.0
    oe8.T_REFLECTION = 180.0
    oe8.T_SOURCE = 0.0
    
    # Elliptical mirror
    oe9.ALPHA = 180.0
    oe9.DUMMY = 0.1
    oe9.FCYL = 1
    oe9.FHIT_C = 1
    oe9.FMIRR = 2
    oe9.F_DEFAULT = 0
    oe9.RLEN1 = 150.0
    oe9.RLEN2 = 150.0
    oe9.RWIDX1 = 15.0
    oe9.RWIDX2 = 15.0
    oe9.SIMAG = 1700.0
    oe9.SSOUR = 8800.0
    oe9.THETA = 88.9056506113
    oe9.T_IMAGE = 175.0
    oe9.T_INCIDENCE = 88.9056506113
    oe9.T_REFLECTION = 88.9056506113
    oe9.T_SOURCE = 8800.0
    
    # Cylindrical mirror
    oe10.ALPHA = 180.0
    oe10.CIL_ANG = 90.0
    oe10.DUMMY = 0.1
    oe10.FCYL = 1
    oe10.FHIT_C = 1
    oe10.FMIRR = 1
    oe10.FWRITE = 1
    oe10.F_EXT = 1
    oe10.RLEN1 = 100.0
    oe10.RLEN2 = 100.0
    oe10.RMIRR = 51.8
    oe10.RWIDX1 = 5.0
    oe10.RWIDX2 = 5.0
    oe10.T_IMAGE = 1350.0
    oe10.T_INCIDENCE  = 90.0 - theta_mrad * 1e-3 * 180 / numpy.pi # 88.9056506113
    oe10.T_REFLECTION = 90.0 - theta_mrad * 1e-3 * 180 / numpy.pi # 88.9056506113
    oe10.T_SOURCE = 175.0
    
    
    
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
    
    
    #
    #run optical element 2
    #
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")
    
    beam.traceOE(oe2,2)
    
    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")
    
    
    #
    #run optical element 3
    #
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")
    
    beam.traceOE(oe3,3)
    
    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")
    
    
    #
    #run optical element 4
    #
    print("    Running optical element: %d"%(4))
    if iwrite:
        oe4.write("start.04")
    
    beam.traceOE(oe4,4)
    
    if iwrite:
        oe4.write("end.04")
        beam.write("star.04")
    
    
    #
    #run optical element 5
    #
    print("    Running optical element: %d"%(5))
    if iwrite:
        oe5.write("start.05")
    
    beam.traceOE(oe5,5)
    
    if iwrite:
        oe5.write("end.05")
        beam.write("star.05")
    
    
    #
    #run optical element 6
    #
    print("    Running optical element: %d"%(6))
    if iwrite:
        oe6.write("start.06")
    
    beam.traceOE(oe6,6)
    
    if iwrite:
        oe6.write("end.06")
        beam.write("star.06")
    
    
    #
    #run optical element 7
    #
    print("    Running optical element: %d"%(7))
    if iwrite:
        oe7.write("start.07")
    
    beam.traceOE(oe7,7)
    
    if iwrite:
        oe7.write("end.07")
        beam.write("star.07")
    
    
    #
    #run optical element 8
    #
    print("    Running optical element: %d"%(8))
    if iwrite:
        oe8.write("start.08")
    
    beam.traceOE(oe8,8)
    
    if iwrite:
        oe8.write("end.08")
        beam.write("star.08")
    
    
    #
    #run optical element 9
    #
    print("    Running optical element: %d"%(9))
    if iwrite:
        oe9.write("start.09")
    
    beam.traceOE(oe9,9)
    
    if iwrite:
        oe9.write("end.09")
        beam.write("star.09")
    
    
    #
    #run optical element 10
    #
    print("    Running optical element: %d"%(10))
    if iwrite:
        oe10.write("start.10")
    
    beam.traceOE(oe10,10)
    
    if iwrite:
        oe10.write("end.10")
        beam.write("star.10")
    
    
    # ~ Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    return beam
    




fnames = [
        '', 
        b'/home/srio/Oasys/ID32/EBS_88_46_polar_H_uz.dat',
        b'/home/srio/Oasys/ID32/EBS_88_46_polar_V_uz.dat',
        b'/home/srio/Oasys/ID32/EBS_68_66_polar_H_uz.dat',
        b'/home/srio/Oasys/ID32/EBS_68_66_polar_V_uz.dat'
        ]

angles = [19.1, 18.875, 19.40, 20.2, 19.45]
            
fig, ax = plt.subplots(2, len(fnames), squeeze=False)
for i, fname in enumerate(fnames):
    beam = trace_rixs_branch(fname, theta_mrad=angles[i])
    nbins = 201
    x = numpy.linspace(-.05, .05, nbins)
    z = numpy.linspace(-.002, .01, nbins)
    
    spot = beam.histo2(1, 3, nbins=nbins, nolost=1, xrange=(x.min(), x.max()), 
        yrange=(z.min(), z.max()), ref=23)
    center = (0, 0)
    ticket = Shadow.ShadowTools.focnew(beam, center=center)
    y = numpy.linspace(-30.0, 30.0, 101) #1001)
    focn = 2.35*Shadow.ShadowTools.focnew_scan(ticket["AX"], y)

    # focn2_h = numpy.zeros_like(focn)
    # focn2_v = numpy.zeros_like(focn)
    # for j in range(focn2_h.size):
    #     beami = beam.duplicate()
    #     beami.retrace(y[j])
    #     fact = 3.0
    #     spot_i = beami.histo2(1, 3, nbins=int(nbins*fact), nolost=1, xrange=(fact*x.min(), fact*x.max()),
    #                        yrange=(fact*z.min(), fact*z.max()), ref=23)
    #     print("j, dist, h, v", j, y[j], spot_i["fwhm_h"], spot_i["fwhm_v"])
    #     focn2_h[j] = spot_i["fwhm_h"]
    #     focn2_v[j] = spot_i["fwhm_v"]

    
    ax[0,i].set_title(str(fname).split('/')[-1])
    ax[0,i].imshow(spot['histogram'].T, extent=(x.min(), x.max(), z.min(), z.max()), origin='lower', aspect='auto')
    ax[1,i].set_title("theta = %5.2f mrad" % angles[i])
    ax[1,i].plot(y, 2.35*Shadow.ShadowTools.focnew_scan(ticket["AX"], y))
    ax[1,i].plot(y, 2.35*Shadow.ShadowTools.focnew_scan(ticket["AZ"], y))
    # ax[1,i].plot(y, focn2_h)
    # ax[1,i].plot(y, focn2_v)

    ax[1,i].set_ylim(0, .05)

plt.show()
    

# fnames = [
#         '',
#         b'C:/Users/kkummer/Desktop/_ID32_ray_tracing/EBS_88_46_polar_H_uz_930eV.dat',
#         b'C:/Users/kkummer/Desktop/_ID32_ray_tracing/EBS_88_46_polar_V_uz_930eV.dat',
#         b'C:/Users/kkummer/Desktop/_ID32_ray_tracing/EBS_68_66_polar_H_uz_930eV.dat',
#         b'C:/Users/kkummer/Desktop/_ID32_ray_tracing/EBS_68_66_polar_V_uz_930eV.dat'
#         ]
#
# fig, ax = plt.subplots(2, len(fnames), squeeze=False)
# for i, fname in enumerate(fnames):
#     beam = trace_rixs_branch(fname)
#     nbins = 201
#     x = numpy.linspace(-.05, .05, nbins)
#     z = numpy.linspace(-.002, .01, nbins)
#
#     spot = beam.histo2(1, 3, nbins=nbins, nolost=1, xrange=(x.min(), x.max()),
#         yrange=(z.min(), z.max()))
#     center = (0, 0)
#     ticket = Shadow.ShadowTools.focnew(beam, center=center)
#     y = numpy.linspace(-20.0, 20.0, 1001)
#     focn = 2.35*Shadow.ShadowTools.focnew_scan(ticket["AX"], y)
#
#     ax[0,i].set_title(str(fname).split('/')[-1])
#     ax[0,i].imshow(spot['histogram'].T, extent=(x.min(), x.max(), z.min(), z.max()), origin='lower', aspect='auto')
#     ax[1,i].plot(y, 2.35*Shadow.ShadowTools.focnew_scan(ticket["AX"], y))
#     axt = plt.twinx(ax[1,i])
#     axt.plot(y, 2.35*Shadow.ShadowTools.focnew_scan(ticket["AZ"], y), color='C3')
#     ax[1,i].set_ylim(0, .05)
#     axt.set_ylim(0, .005)
    
# plt.show()
