from srxraylib.plot.gol import plot
import numpy
si = numpy.loadtxt("rothSi.dat")
ni = numpy.loadtxt("rothNi.dat")
pt = numpy.loadtxt("rothPt.dat")

plot(si[:,0,], si[:,-2],
    ni[:,0,], ni[:,-2],
    pt[:,0,], pt[:,-2],
    si[:, 0, ], si[:, -1],
    legend=["Si","Ni","Pt","Source"])


