import os, numpy

# os.system("rm tmp.inp tmp.out tmp.ene")
ENERGIES = numpy.linspace(1000.0,30000,100)
#
#
# for energy in ENERGIES:
#
#     #energy=10000.0
#     f = open("tmp.inp",'w')
#     f.write("10000\n0.03\n10\n90\n100\n1\n%g" % energy)
#     f.close()
#
#     os.system("python pescao.py < tmp.inp")
#     os.system("cat OutputValues.dat")
#     os.system('grep "Absorbed power fraction in the silicon"  OutputValues.dat >> tmp.out')
#     os.system("echo %d >> tmp.ene" % energy)



rho = 1.848
l= ENERGIES.size

# xraylib
import xraylib
XRL_MU = numpy.zeros(l)
XRL_MU_E = numpy.zeros(l)

for i in range(l):
    XRL_MU[i] =   rho * xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), 1e-3*ENERGIES[i])
    XRL_MU_E[i] = rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), 1e-3*ENERGIES[i])
    # print("     >>",ENERGIES[i],rho*muoverrhos[i,1], XRL_MU[i], XRL_MU_E[i], XRL_MU_E[i] / XRL_MU[i])

from srxraylib.plot.gol import plot
# plot(
#      ENERGIES, XRL_MU,
#      ENERGIES, XRL_MU_E,
#      xlog=1,ylog=1,title="Be",xrange=[1e3,1e5],xtitle="Photon energy [eV]",ytitle="mu [cm-1]",
#      legend=["mu xraylib","energy loss mu"])

thickness = 800e-6
XRL_ABS   = 1.0 - numpy.exp(-XRL_MU * thickness * 100)
XRL_ABS_E = 1.0 - numpy.exp(-XRL_MU_E * thickness * 100)

a_ene = numpy.loadtxt("tmp_0p8.ene")
a_out = numpy.loadtxt("tmp_0p8.out")
plot(
     ENERGIES, XRL_ABS,
     ENERGIES, XRL_ABS_E,
     a_ene, 1e-2 * a_out[:,0],
     xlog=0,ylog=0,xrange=[0,20000],title="Be thickness=%g um" %(1e6*thickness),xtitle="Photon energy [eV]",ytitle="Absorption fraction",
     legend=["mu xraylib","energy loss mu",'Monte Carlo'],show=1)

c8 = (1e-2 * a_out[:,0] - XRL_ABS) / (XRL_ABS_E - XRL_ABS)

thickness = 300e-6
XRL_ABS   = 1.0 - numpy.exp(-XRL_MU * thickness * 100)
XRL_ABS_E = 1.0 - numpy.exp(-XRL_MU_E * thickness * 100)

print("Be")
energy = 5.0
print("Density: %g g/cm3" % rho)
print("Be thickness: %g m" % thickness)
print("Energy: %g keV"      % energy)
print("mu_a/ro: %g cm2/g"     %  xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), energy))
print("mu_en/ro: %g cm2/q"    %  xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), energy))
print("mu_a/ro: %g cm-1"     %  (rho *  xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), energy)))
print("mu_en/ro: %g cm-1"    %  (rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), energy)))
print("1-exp(-mu*t): %g "     %  (1.0 - numpy.exp(- thickness * 100 * rho *  xraylib.CS_Total(xraylib.SymbolToAtomicNumber("Be"), energy))))
print("1-exp(-mu*t): %g "     %  (1.0 - numpy.exp(- thickness * 100 * rho * xraylib.CS_Energy(xraylib.SymbolToAtomicNumber("Be"), energy))))
# a_ene = numpy.loadtxt("tmp_0p3.ene")
# a_out = numpy.loadtxt("tmp_0p3.out")
# plot(
#      ENERGIES, XRL_ABS,
#      ENERGIES, XRL_ABS_E,
#      a_ene, 1e-2 * a_out[:,0],
#      xlog=0,ylog=0,xrange=[0,20000],title="Be thickness=%g um" %(1e6*thickness),xtitle="Photon energy [eV]",ytitle="Absorption fraction",
#      legend=["mu xraylib","energy loss mu",'Monte Carlo'], show=0)
#
# c3 = (1e-2 * a_out[:,0] - XRL_ABS) / (XRL_ABS_E - XRL_ABS)
#
#
#
#
# thickness = 100e-6
# XRL_ABS   = 1.0 - numpy.exp(-XRL_MU * thickness * 100)
# XRL_ABS_E = 1.0 - numpy.exp(-XRL_MU_E * thickness * 100)
# a_ene = numpy.loadtxt("tmp_0p1.ene")
# a_out = numpy.loadtxt("tmp_0p1.out")
# plot(
#      ENERGIES, XRL_ABS,
#      ENERGIES, XRL_ABS_E,
#      a_ene, 1e-2 * a_out[:,0],
#      xlog=0,ylog=0,xrange=[0,20000],title="Be thickness=%g um" %(1e6*thickness),xtitle="Photon energy [eV]",ytitle="Absorption fraction",
#      legend=["mu xraylib","energy loss mu",'Monte Carlo'], show=0)
#
# c3 = (1e-2 * a_out[:,0] - XRL_ABS) / (XRL_ABS_E - XRL_ABS)
#
# plot(ENERGIES,c3, ENERGIES, c8,
#      legend=["c3","c8"], ylog=1)