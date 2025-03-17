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
    try:
        XRL_MU[i] = rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), muoverrhos[i,0])
        XRL_MU_E[i] = rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), muoverrhos[i,0])
    except:
        XRL_MU[i] = 0
        XRL_MU_E[i] = 0
        # print("     >>",muoverrhos[i,0],rho*muoverrhos[i,1], XRL_MU[i], XRL_MU_E[i], XRL_MU_E[i] / XRL_MU[i])

from srxraylib.plot.gol import plot
# plot(muoverrhos[:,0], rho*muoverrhos[:,1],
#      muoverrhos[:,0], XRL_MU,
#      muoverrhos[:, 0], XRL_MU_E,
#      xlog=1,ylog=1,title="Be",xrange=[1,1e3],xtitle="Photon energy [keV]",ytitle="mu [cm-1]",
#      legend=["mu xro","mu xraylib","energy loss mu"])





nist = [
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
    [2.00000E+01, 1.227E-02, 9.294E-03]]
nist = numpy.array(nist)
print(">>", nist.shape)

plot(nist[:,0], nist[:,1],
     nist[:,0], nist[:,2],
     xlog=1,ylog=1,title="Be",xtitle="Photon energy [MeV]",ytitle="mu [cm2/g]",
     legend=["mu nist", "mu_en nist"])

plot(muoverrhos[:,0], XRL_MU/rho,
     muoverrhos[:, 0], XRL_MU_E/rho,
     nist[:,0]*1e+3, nist[:,1],
     nist[:,0]*1e+3, nist[:,2],
     xlog=1,ylog=1,title="Be",xrange=[1,1e3],xtitle="Photon energy [MeV]",ytitle="mu [cm2/g]",
     legend=["mu xraylib","mu_en xraylib", "mu nist", "mu_en nist"])

