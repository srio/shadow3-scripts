import numpy


# from rxoptics.de
rho=1.85
muoverrhos= [[1.0,604.1],[1.5,179.7],[2.0,74.69],[3.0,21.27],[4.0,8.685],[5.0,4.369],[6.0,2.527],[8.0,1.124],[10.0,0.6466],[15.0,0.307],[20.0,0.2251],[30.0,0.1792],[40.0,0.164],[50.0,0.1554],[60.0,0.1493],[80.0,0.1401],[100.0,0.1328],[150.0,0.119],[200.0,0.1089],[300.0,0.09463],[400.0,0.08471],[500.0,0.07739],[600.0,0.07155],[800.0,0.06286],[1000.0,0.05652],[1250.0,0.05054],[1500.0,0.04597],[2000.0,0.03938],[3000.0,0.03138],[4000.0,0.02664],[5000.0,0.02347],[6000.0,0.02121],[8000.0,0.01819],[10000.0,0.01627],[15000.0,0.01361],[20000.0,0.01227]]
muoverrhos = numpy.array(muoverrhos)
l= muoverrhos.shape[0]

# xraylib
import xraylib
XRL_MU = numpy.zeros(l)
XRL_MU_E = numpy.zeros(l)

for i in range(l):
    XRL_MU[i] = rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), muoverrhos[i,0])
    XRL_MU_E[i] = rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), muoverrhos[i,0])
    print("     >>",muoverrhos[i,0],rho*muoverrhos[i,1], XRL_MU[i], XRL_MU_E[i], XRL_MU_E[i] / XRL_MU[i])

from srxraylib.plot.gol import plot
plot(muoverrhos[:,0], rho*muoverrhos[:,1],
     muoverrhos[:,0], XRL_MU,
     muoverrhos[:, 0], XRL_MU_E,
     xlog=1,ylog=1,title="Be",xrange=[1,1e3],xtitle="Photon energy [keV]",ytitle="mu [cm-1]",
     legend=["mu xro","mu xraylib","energy loss mu"])

