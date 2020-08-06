"""
File with trajectory written to file: /users/srio/Oasys/tmp.traj

wiggler_cdf: Electron beam energy (from velocities) = 3.000355 GeV

wiggler_cdf: gamma (from velocities) = 5870.853556 GeV
wiggler_cdf: Curvature (min)) = 0.000000 m^-1
wiggler_cdf:           (max)    0.199920 m^-1
wiggler_cdf: Radius of curvature (max) = 81689012171814624.000000 m
wiggler_cdf:                     (min) = 5.002009 m
wiggler_cdf: Critical Energy (max.) = 11973.937061 eV
wiggler_cdf:                 (min.) = 0.000000 eV
wiggler_cdf: Total no.of photons = 1.690471e+17 (in DE=99900.000 eV)
wiggler_cdf: File with wiggler cdf written to file: b'/users/srio/Oasys/xshwig.sha'

Electron beam energy (from velocities) = 3.000355 GeV

gamma (from velocities) = 5870.851896
curvature (max) = 0.199920 m
          (min) = 0.000000 m
Radius of curvature (max) = 81689012171830928.000000 m
                    (min) = 5.002009 m
Critical Energy (max.) = 11973.926903 eV
                (min.) = 0.000000 eV
File with wiggler spectrum written to file: spectrum.dat

Total power (from integral of spectrum): 10106.973910 W

Total number of photons (from integral of spectrum): 1.62115e+19

"""

#
# script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
#
from srxraylib.sources import srfunc
from srxraylib.plot.gol import plot, plot_image, plot_scatter
import numpy

(traj, pars) = srfunc.wiggler_trajectory(
    b_from = 0,
    inData = "",
    nPer = 3, #37,
    nTrajPoints = 501,
    ener_gev = 3.0,
    per = 0.12,
    kValue = 22.416,
    trajFile = "tmp.traj",
    shift_x_flag = 0,
    shift_x_value = 0.0,
    shift_betax_flag = 0,
    shift_betax_value = 0.0)



#
# calculate cdf and write file for Shadow/Source
#

# srfunc.wiggler_cdf(traj,
#     enerMin = 100.0,
#     enerMax = 100000.0,
#     enerPoints = 1001,
#     outFile = b'/users/srio/Oasys/xshwig.sha',
#     elliptical = False)



# propagate horizontal divergence

distance = 30.0

Y = traj[1,:].copy()
divX = traj[3,:].copy()
curX = traj[6,:].copy()
posX = divX * (distance + Y)

# detX = numpy.linspace(posX.min(), posX.max(), divX.size)
plot(Y,numpy.abs(curX)**2, xtitle="Y", ytitle="curX")
plot(divX, curX**2, xtitle="divX", ytitle="Intensity")
plot_scatter(posX, curX**2, xtitle="detector X", ytitle="intensity")

# hisX, edgesX = numpy.histogram(posX, bins=100, range=None, weights=curX**2, density=False)
# hisX0, edgesX = numpy.histogram(posX, bins=100, range=None, weights=None, density=False)
# plot(numpy.arange(hisX.size), hisX/hisX0)

# enerMin = 100.0,
# enerMax = 100000.0,
# nPoints = 500,
#
# e, f, w = srfunc.wiggler_spectrum(traj,
#     enerMin = enerMin,
#     enerMax = enerMax,
#     nPoints = nPoints,
#     electronCurrent = 100.0 * 1e-3,
#     outFile = "spectrum.dat",
#     elliptical = False)
#
#
#
# # plot(e, f, xlog=False, ylog=False, show=False,
# # xtitle = "Photon energy [eV]", ytitle = "Flux [Photons/s/0.1%bw]", title = "Flux")
# # plot(e, w, xlog=False, ylog=False, show=True,
# # xtitle = "Photon energy [eV]", ytitle = "Spectral Power [E/eV]", title = "Spectral Power")
#
#
#
# # :return:           (traj,pars)
# #        traj: a variable to store the output matrix (8 colums with:
# #        x[m]  y[m]  z[m]  BetaX  BetaY  BetaZ  Curvature  B[T] )
# #        pars: a variable with text info
#
# # plot(traj[1,:], traj[0,:], xtitle="Y", ytitle="X")
# # plot(traj[1,:], traj[3,:], xtitle="Y", ytitle="betaX")
# # print(traj.shape)
