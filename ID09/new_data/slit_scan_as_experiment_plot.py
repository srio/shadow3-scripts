import numpy
from srxraylib.plot.gol import plot



a = numpy.loadtxt("slit_scan_as_experiment.dat")

err = numpy.random.normal(loc=0.0, scale=2.8e-6/2.355, size=a.shape[0]) * 1
# print(err)
a[:,1] += err

afit = numpy.polyfit(a[:,0], a[:,1], 1)
print(afit)

a2 = a[:,0] * afit[0] + afit[1]

plot(1e3 * a[:,0], a[:,1],
     1e3 * a[:,0], a2,
     1e3 * a[:,0], a[:,1] - a2,
     xtitle="Slit position [mm]",
     ytitle="Position on xeye [m]",
     legend=["with slope error", "linear fit", "difference"],
     show=1
    )

# slope = (a1[:,1] - a0[:,1]) / (55.2290 - 54.359) * (-1)

slope = (a[:,1] - a2) / (55.2290 - 54.359) * (-1) / 2

print("Slope error RMS: ", numpy.std(slope))

plot(1e3 * a[:,0], slope,
     xtitle="Slit position [mm]",
     ytitle="slope on mirror",
     title="retrieved slope error",
     show=1)


plot(1e3 * a[:,0], numpy.cumsum(slope) ,
     xtitle="Slit position [mm]",
     ytitle="height on mirror",
     title="retrieved height error",
     show=1)


#
# a0 = numpy.loadtxt("slit_scan_as_experiment0.dat")
# a1 = numpy.loadtxt("slit_scan_as_experiment1.dat")
#
# err = numpy.random.normal(loc=0.0, scale=2e-6/2.355, size=a1.shape[0]) * 1
# # print(err)
# a1[:,1] += err
#
# afit = numpy.polyfit(a1[:,0], a1[:,1], 1)
# print(afit)
# a2 = a1[:,0] * afit[0] + afit[1]
#
# plot(1e3 * a0[:,0], a0[:,1],
#      1e3 * a1[:,0], a1[:,1],
#      1e3 * a1[:,0], a2,
#      1e3 * a0[:,0], a1[:,1] - a0[:,1],
#      xtitle="Slit position [mm]",
#      ytitle="Position on xeye [m]",
#      legend=["no slope error","with slope error", "linear fit", "difference"],
#      show=1
#     )
#
# # slope = (a1[:,1] - a0[:,1]) / (55.2290 - 54.359) * (-1)
#
# slope = (a1[:,1] - a2) / (55.2290 - 54.359) * (-1)
#
#
#
# plot(1e3 * a0[:,0], slope,
#      xtitle="Slit position [mm]",
#      ytitle="slope on mirror",
#      title="retrieved slope error",
#      show=1)
#
# plot(1e3 * a0[:,0], numpy.cumsum(slope) ,
#      xtitle="Slit position [mm]",
#      ytitle="height on mirror",
#      title="retrieved height error",
#      show=1)

# coordinate_on_mirror = numpy.linspace(a[:,1].min(), a[:,1].max(), a.shape[0])
#
# plot(1e3 * coordinate_on_mirror, slope,
#      xtitle="Mirror position [m]",
#      ytitle="slope on mirror",
#      title="retrieved slope error",
#      show=1)

# plot(a[:,0], numpy.cumsum(slope, a[:,0]), title="heights")