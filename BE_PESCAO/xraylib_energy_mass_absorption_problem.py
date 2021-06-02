import numpy
import xraylib

from srxraylib.plot.gol import plot, set_qt
set_qt()



#
# Be
#

# energy Mev, mu/rho, mu_en/rho
nist_Be = numpy.array([
    [1.00000E-03, 6.041E+02, 6.035E+02],
    [1.50000E-03, 1.797E+02, 1.791E+02],
    [2.00000E-03, 7.469E+01, 7.422E+01],
    [3.00000E-03, 2.127E+01, 2.090E+01],
    [4.00000E-03, 8.685E+00, 8.367E+00],
    [5.00000E-03, 4.369E+00, 4.081E+00],
    [6.00000E-03, 2.527E+00, 2.260E+00],
    [8.00000E-03, 1.124E+00, 8.839E-01],
    [1.00000E-02, 6.466E-01, 4.255E-01],
    [1.50000E-02, 3.070E-01, 1.143E-01],
    [2.00000E-02, 2.251E-01, 4.780E-02],
    [3.00000E-02, 1.792E-01, 1.898E-02],
    [4.00000E-02, 1.640E-01, 1.438E-02],
    [5.00000E-02, 1.554E-01, 1.401E-02],
    [6.00000E-02, 1.493E-01, 1.468E-02],
    [8.00000E-02, 1.401E-01, 1.658E-02],
    [1.00000E-01, 1.328E-01, 1.836E-02],
    [1.50000E-01, 1.190E-01, 2.157E-02],
    [2.00000E-01, 1.089E-01, 2.353E-02],
    [3.00000E-01, 9.463E-02, 2.548E-02],
    [4.00000E-01, 8.471E-02, 2.620E-02],
    [5.00000E-01, 7.739E-02, 2.639E-02],
    [6.00000E-01, 7.155E-02, 2.627E-02],
    [8.00000E-01, 6.286E-02, 2.565E-02],
    [1.00000E+00, 5.652E-02, 2.483E-02],
    [1.25000E+00, 5.054E-02, 2.373E-02],
    [1.50000E+00, 4.597E-02, 2.268E-02],
    [2.00000E+00, 3.938E-02, 2.083E-02],
    [3.00000E+00, 3.138E-02, 1.806E-02],
    [4.00000E+00, 2.664E-02, 1.617E-02],
    [5.00000E+00, 2.347E-02, 1.479E-02],
    [6.00000E+00, 2.121E-02, 1.377E-02],
    [8.00000E+00, 1.819E-02, 1.233E-02],
    [1.00000E+01, 1.627E-02, 1.138E-02],
    [1.50000E+01, 1.361E-02, 1.001E-02],
    [2.00000E+01, 1.227E-02, 9.294E-03]])


energy_ev = nist_Be[:,0] * 1e6
# xraylib
XRL_MU = numpy.zeros_like(energy_ev)
XRL_MU_E = numpy.zeros_like(energy_ev)

for i in range(energy_ev.size):
    XRL_MU[i] = xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), 1e-3 * energy_ev[i])
    XRL_MU_E[i] = xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), 1e-3 * energy_ev[i])

plot(energy_ev, XRL_MU,
     nist_Be[:,0]*1e+6, nist_Be[:,1],
     xlog=1,ylog=1,title="Be",xtitle="Photon energy [eV]",ytitle="mu",xrange=[1e3,1e6],
     legend=["Be mu xraylib", "Be mu nist"])
plot(
     energy_ev, XRL_MU_E,
     nist_Be[:,0]*1e+6, nist_Be[:,2],
     xlog=1,ylog=1,title="Be",xtitle="Photon energy [eV]",ytitle="mu",xrange=[1e3,1e6],
     legend=["Be mu_en xraylib", "Be mu_en nist"])


nist_C = numpy.array([
    [1.00000E-03, 2.211E+03, 2.209E+03],
    [1.50000E-03, 7.002E+02, 6.990E+02],
    [2.00000E-03, 3.026E+02, 3.016E+02],
    [3.00000E-03, 9.033E+01, 8.963E+01],
    [4.00000E-03, 3.778E+01, 3.723E+01],
    [5.00000E-03, 1.912E+01, 1.866E+01],
    [6.00000E-03, 1.095E+01, 1.054E+01],
    [8.00000E-03, 4.576E+00, 4.242E+00],
    [1.00000E-02, 2.373E+00, 2.078E+00],
    [1.50000E-02, 8.071E-01, 5.627E-01],
    [2.00000E-02, 4.420E-01, 2.238E-01],
    [3.00000E-02, 2.562E-01, 6.614E-02],
    [4.00000E-02, 2.076E-01, 3.343E-02],
    [5.00000E-02, 1.871E-01, 2.397E-02],
    [6.00000E-02, 1.753E-01, 2.098E-02],
    [8.00000E-02, 1.610E-01, 2.037E-02],
    [1.00000E-01, 1.514E-01, 2.147E-02],
    [1.50000E-01, 1.347E-01, 2.449E-02],
    [2.00000E-01, 1.229E-01, 2.655E-02],
    [3.00000E-01, 1.066E-01, 2.870E-02],
    [4.00000E-01, 9.546E-02, 2.950E-02],
    [5.00000E-01, 8.715E-02, 2.969E-02],
    [6.00000E-01, 8.058E-02, 2.956E-02],
    [8.00000E-01, 7.076E-02, 2.885E-02],
    [1.00000E+00, 6.361E-02, 2.792E-02],
    [1.25000E+00, 5.690E-02, 2.669E-02],
    [1.50000E+00, 5.179E-02, 2.551E-02],
    [2.00000E+00, 4.442E-02, 2.345E-02],
    [3.00000E+00, 3.562E-02, 2.048E-02],
    [4.00000E+00, 3.047E-02, 1.849E-02],
    [5.00000E+00, 2.708E-02, 1.710E-02],
    [6.00000E+00, 2.469E-02, 1.607E-02],
    [8.00000E+00, 2.154E-02, 1.468E-02],
    [1.00000E+01, 1.959E-02, 1.380E-02],
    [1.50000E+01, 1.698E-02, 1.258E-02],
    [2.00000E+01, 1.575E-02, 1.198E-02],
 ])

energy_ev = nist_C[:,0] * 1e6
# xraylib
XRL_MU = numpy.zeros_like(energy_ev)
XRL_MU_E = numpy.zeros_like(energy_ev)

for i in range(energy_ev.size):
    XRL_MU[i] = xraylib.CS_Total(xraylib.SymbolToAtomicNumber("C"), 1e-3 * energy_ev[i])
    XRL_MU_E[i] = xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("C"), 1e-3 * energy_ev[i])

plot(energy_ev, XRL_MU,
     nist_C[:,0]*1e+6, nist_C[:,1],
     xlog=1,ylog=1,title="C",xtitle="Photon energy [eV]",ytitle="mu",xrange=[1e3,1e6],
     legend=["C mu xraylib", "C mu nist"])
plot(
     energy_ev, XRL_MU_E,
     nist_C[:,0]*1e+6, nist_C[:,2],
     xlog=1,ylog=1,title="C",xtitle="Photon energy [eV]",ytitle="mu",xrange=[1e3,1e6],
     legend=["C mu_en xraylib", "C mu_en nist"])