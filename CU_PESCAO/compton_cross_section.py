import numpy
import xraylib

from srxraylib.plot.gol import plot, set_qt
set_qt()





energy_kev = numpy.linspace(1,1000,1000)
# xraylib
Cu_MU = numpy.zeros_like(energy_kev)
Cu_MU_C = numpy.zeros_like(energy_kev)
Be_MU = numpy.zeros_like(energy_kev)
Be_MU_C = numpy.zeros_like(energy_kev)

for i in range(energy_kev.size):
    Cu_MU[i] =   xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Si"), energy_kev[i])
    Cu_MU_C[i] = xraylib.CS_Compt(xraylib.SymbolToAtomicNumber("Si"), energy_kev[i])
    Be_MU[i] =   xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), energy_kev[i])
    Be_MU_C[i] = xraylib.CS_Compt(xraylib.SymbolToAtomicNumber("Be"), energy_kev[i])

plot(   energy_kev, Cu_MU,
        energy_kev, Cu_MU_C,
        energy_kev, Be_MU,
        energy_kev, Be_MU_C,
        xlog=1,ylog=1,title="Cu",xtitle="Photon energy [keV]",ytitle="mu",
        legend=["Si mu",
                "Si mu (Compton)",
                "Be mu",
                "Be mu (Compton)"],
        linestyle=[':', None, ':', None],
        color=['r', 'r', 'b', 'b'],
        show=0)

print("Integrated Cu, Be: ", Cu_MU_C.sum(), Be_MU_C.sum())


angle = numpy.linspace(0, 2*numpy.pi, 500)
Si_ang_Compton = numpy.zeros_like(angle)
Si_ang_Rayleigh = numpy.zeros_like(angle)
Be_ang_Compton = numpy.zeros_like(angle)
Be_ang_Rayleigh = numpy.zeros_like(angle)
OO_ang_Compton = numpy.zeros_like(angle)
OO_ang_Rayleigh = numpy.zeros_like(angle)



energy = 50.0
for i in range(angle.size):
    Si_ang_Compton[i] =  xraylib.DCSb_Compt(xraylib.SymbolToAtomicNumber("Si"), energy, angle[i])
    Si_ang_Rayleigh[i] = xraylib.DCSb_Rayl (xraylib.SymbolToAtomicNumber("Si"), energy, angle[i])
    Be_ang_Compton[i] =  xraylib.DCSb_Compt(xraylib.SymbolToAtomicNumber("Be"), energy, angle[i])
    Be_ang_Rayleigh[i] =  xraylib.DCSb_Rayl(xraylib.SymbolToAtomicNumber("Be"), energy, angle[i])
    OO_ang_Compton[i] =  xraylib.DCS_KN  (energy, angle[i])
    OO_ang_Rayleigh[i] = xraylib.DCS_Thoms(       angle[i])

print("Integrated Angle Si, Be: ", Si_ang_Compton.sum(), Be_ang_Compton.sum())

plot(angle, Si_ang_Compton,
     angle, Si_ang_Rayleigh,
     angle, Be_ang_Compton,
     angle, Be_ang_Rayleigh,
     legend=["Compton  Si",
             "Rayleigh Si",
             "Compton  Be",
             "Rayleigh Be"],
     linestyle=[None, ':', None, ':'],
     color=['r', 'r', 'b', 'b'],
     show=0)

import numpy as np
import matplotlib.pyplot as plt


# r = np.arange(0, 2, 0.01)
# theta = angle # 2 * np.pi * r
# r = Cu_ang_Compton

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(angle, Si_ang_Compton,  color='r', label="Si")
ax.plot(angle, Be_ang_Compton,  color='b', label="Be")
ax.plot(angle, OO_ang_Compton,  color='g', label="1e")

ax.legend()
# ax.plot(angle, OO_ang_Rayleigh, color='k', label="Cu")
# ax.set_rmax(2)
ax.set_rticks(numpy.linspace(0,(numpy.concatenate((Si_ang_Compton, Be_ang_Compton))).max(),10))  # Less radial ticks
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("Energy=%g keV" % energy, va='bottom')
plt.show()

