import numpy

from srxraylib.plot.gol import plot, set_qt

set_qt()

a = numpy.loadtxt("Dia")

import xraylib

print(">>>>", a.shape)
b = numpy.zeros((a.shape[0], 3))

for i in range(a.shape[0]):
    b[i,0] = a[i,0]
    n = xraylib.Refractive_Index("C", a[i,0] * 1e-3, 3.51)
    b[i,1] = 1.0 - n.real
    b[i,2] = n.imag



plot(a[:, 0], a[:, 1],
     a[:, 0], a[:, 2],
     b[:, 0], b[:, 1],
     b[:, 0], b[:, 2],
     legend=["srcalc delta",
             "srcalc beta",
             "xraylib delta",
             "xraylib beta"],
     xlog=1,ylog=1,
     xtitle='Photon energy [eV]', ytitle="delta or beta")


f = open("Dia-xraylib",'w')
for i in range(a.shape[0]):
    f.write("%g %g %g\n" % (b[i,0], b[i,1], b[i,2]) )

f.close()
print("File written to disk: Dia-xraylib")