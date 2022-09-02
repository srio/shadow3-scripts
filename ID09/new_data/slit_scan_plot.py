import numpy
from srxraylib.plot.gol import plot

a = numpy.loadtxt("slit_scan.dat")
plot(1e3 * a[:,0], a[:,1],
     1e3 * a[:,0], a[:,2],
     1e3 * a[:,0], a[:,2] - a[:,1],
     xtitle="Slit position [mm]",
     ytitle="Position on mirror surface [m]",
     legend=["no slope error","with slope error", "difference"],
     show=0
    )

slope = (a[:,2] - a[:,1]) * 2.5e-3 / (55.2290  - 44.5400 ) / 2

print("Slope error RMS: ", numpy.std(slope))

plot(1e3 * a[:,0], slope,
     xtitle="Slit position [mm]",
     ytitle="slope on mirror",
     title="retrieved slope error",
     show=1)

coordinate_on_mirror = numpy.linspace(a[:,1].min(), a[:,1].max(), a.shape[0])

plot(1e3 * coordinate_on_mirror, slope,
     xtitle="Mirror position [m]",
     ytitle="slope on mirror",
     title="retrieved slope error",
     show=1)

plot(1e3 * coordinate_on_mirror, numpy.cumsum(slope) * (coordinate_on_mirror[1]-coordinate_on_mirror[0]),
     xtitle="Mirror position [m]",
     ytitle="heights on mirror",
     title="retrieved slope error",
     show=1)

# plot(a[:,0], numpy.cumsum(slope, a[:,0]), title="heights")