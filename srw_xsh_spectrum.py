#############################################################################
# SRWLIB Example#6: Calculating spectral flux of undulator radiation 
# by finite-emittance electron beam collected through a finite aperture
# and power density distribution of this radiation (integrated over all photon energies)
# v 0.02
#
# Modified by srio@esrf.eu 2012-12-19
# set the inputs easier and dump all results (flux, power density and trajectory 
# in a single spec-formatted file
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
#import os
#import sys

#**********************Predefined Input Parameters:




import json
inFileTxt = "xshundul.json"
with open(inFileTxt, mode='r') as f1:
    h = json.load(f1)

#
# list all non-empty keywords
#
#print "-----------------------------------------------------"
for i,j in h.items():
    if (j != None):
        print ("%s = %s" % (i,j) )
#print "-----------------------------------------------------"
#print "k: ",h['K']


# as defined in glossary (Gaussian Beam)
ElectronEnergy = h['E_ENERGY']
ElectronCurrent = h['INTENSITY']
ElectronBeamSizeH = h['SX']*1e-2
ElectronBeamSizeV = h['SZ']*1e-2
if (ElectronBeamSizeH > 0):
    ElectronBeamDivergenceH = h['EX']*1e-2/ElectronBeamSizeH
else:
    ElectronBeamDivergenceH = 0e0

if (ElectronBeamSizeV > 0):
    ElectronBeamDivergenceV = h['EZ']*1e-2/ElectronBeamSizeV
else:
    ElectronBeamDivergenceV = 0e0

print (">> ElectronEnergy [GeV]: %f"%(ElectronEnergy) )
print (">> ElectronCurrent [A]: %f"%(ElectronCurrent) )
print (">> ElectronBeamSizeH [m]: %f"%(ElectronBeamSizeH) )
print (">> ElectronBeamSizeV [m]: %f"%(ElectronBeamSizeV) )
print (">> ElectronBeamDivergenceH [rad]: %f"%(ElectronBeamDivergenceH) )
print (">> ElectronBeamDivergenceV [rad]: %f"%(ElectronBeamDivergenceV) )

# as defined in glossary (Insertion Device)
PeriodID = h['LAMBDAU']
N = h['NPERIODS']
print (">> PeriodID [m]: %f"%(PeriodID) )
print (">> N: %d"%(N) )
#case 1 E_III = 15816
#Kv = 1.452
#outFile = 'srw_id28_case1.spec'
#case 3 E_III = 21747
#Kv = 0.994
#outFile = 'srw_id28_case3.spec'
#case 4 E_V = 21747
#Kv = 1.7263
#outFile = 'srw_id28_case4.spec'
#case 5 E_V = 23750
Kv = h['K'] # 1.6
print (">> Kv: %f"%(Kv) )
outFile = 'srw_xsh_spectrum.spec'

#case 2 E_II = 15816
#Kv = 0.46
#PeriodID = 17.6e-3
#N = 91
#outFile = 'srw_id28_case2.spec'


# as defined in glossary (driftspace) = distance from undulator to slit
d = 50.0
# as defined in glossary (slit)
gapH = d * h['MAXANGLE']*1e-3 * 2.0
gapV = gapH
#others (scan, output)
PhotonEnergyMin = h['EMIN']
PhotonEnergyMax = h['EMAX']
print (">> PhotonEnergyMin [eV]: %f"%(PhotonEnergyMin) )
print (">> PhotonEnergyMax [eV]: %f"%(PhotonEnergyMax) )
print (">> slit gapH [um]: %f"%(gapH*1e6) )
print (">> slit gapV [um]: %f"%(gapV*1e6) )
print (">> slit distance d [m]: %f"%(d) )
#PhotonEnergyPoints = 1 + int( (PhotonEnergyMax-PhotonEnergyMin)/h['ESTEP'] ) #10 #0
PhotonEnergyPoints = h['NPOINTS_ENERGY']
AnglePoints = h['NPOINTS_ANGLE']
print (">> PhotonEnergyPoints [eV]: %d"%(PhotonEnergyPoints) )
Nmax = 10 # 24 # 91

#derived
B0 = Kv/0.934/(PeriodID*1e2)

#**********************End Predefined Input Parameters:


#print('SRWLIB Python Example # 6:')
print('Running SRW (SRWLIB Python)')
print('Calculating spectral flux of undulator radiation by finite-emittance electron beam collected through a finite aperture and power density distribution of this radiation (integrated over all photon energies)')

#***********Undulator
harmB = SRWLMagFldH() #magnetic field harmonic
harmB.n = 1 #harmonic number
harmB.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB.B = B0 #magnetic field amplitude [T]
und = SRWLMagFldU([harmB])
und.per = PeriodID #period length [m]
und.nPer = N #number of periods (will be rounded to integer)
magFldCnt = SRWLMagFldC([und], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all magnetic field elements

#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = ElectronCurrent #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = ElectronEnergy/0.51099890221e-03 #relative energy
sigEperE = 1e-8 # 0.00089 #relative RMS energy spread
sigX =  ElectronBeamSizeH #horizontal RMS size of e-beam [m]
sigXp = ElectronBeamDivergenceH #horizontal RMS angular divergence [rad]
sigY =  ElectronBeamSizeV #vertical RMS size of e-beam [m]
sigYp = ElectronBeamDivergenceV #vertical RMS angular divergence [rad]
#2nd order stat. moments:
eBeam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2> 
eBeam.arStatMom2[1] = 0 #<(x-<x>)(x'-<x'>)>
eBeam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2> 
eBeam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
eBeam.arStatMom2[4] = 0 #<(y-<y>)(y'-<y'>)>
eBeam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
eBeam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E>^2

#***********Precision Parameters
arPrecF = [0]*5 #for spectral flux vs photon energy
arPrecF[0] = 1 #initial UR harmonic to take into account
arPrecF[1] = Nmax #final UR harmonic to take into account
arPrecF[2] = 1.5 #longitudinal integration precision parameter
arPrecF[3] = 1.5 #azimuthal integration precision parameter
arPrecF[4] = 1 #calculate flux (1) or flux per unit surface (2)

arPrecP = [0]*5 #for power density
arPrecP[0] = 1.5 #precision factor
arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

#***********UR Stokes Parameters (mesh) for Spectral Flux
stkF = SRWLStokes() #for spectral flux vs photon energy
#srio stkF.allocate(10000, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
stkF.allocate(PhotonEnergyPoints, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
stkF.mesh.zStart = d #longitudinal position [m] at which UR has to be calculated
stkF.mesh.eStart = PhotonEnergyMin #initial photon energy [eV]
stkF.mesh.eFin =   PhotonEnergyMax #final photon energy [eV]
stkF.mesh.xStart = -gapH/2 #initial horizontal position [m]
stkF.mesh.xFin =    gapH/2 #final horizontal position [m]
stkF.mesh.yStart = -gapV/2 #initial vertical position [m]
stkF.mesh.yFin =    gapV/2 #final vertical position [m]

stkP = SRWLStokes() #for power density
stkP.allocate(1, 51, 51) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
stkP.mesh.zStart = d #longitudinal position [m] at which power density has to be calculated
stkP.mesh.xStart = -gapH/2 #initial horizontal position [m]
stkP.mesh.xFin =    gapH/2 #final horizontal position [m]
stkP.mesh.yStart = -gapV/2 #initial vertical position [m]
stkP.mesh.yFin =    gapV/2 #final vertical position [m]

#sys.exit(0)

#**********************Calculation (SRWLIB function calls)
print('   Performing Spectral Flux (Stokes parameters) calculation ... ', end='')
srwl.CalcStokesUR(stkF, eBeam, und, arPrecF)
print('done')

partTraj = SRWLPrtTrj() #defining auxiliary trajectory structure
partTraj.partInitCond = eBeam.partStatMom1
partTraj.allocate(20001) 
partTraj.ctStart = -1.6 #Start Time for the calculation
partTraj.ctEnd = 1.6
partTraj = srwl.CalcPartTraj(partTraj, magFldCnt, [1])

#print('   Performing Power Density calculation (from trajectory) ... ', end='')
#srwl.CalcPowDenSR(stkP, eBeam, partTraj, 0, arPrecP)
#print('done')

print('   Performing Power Density calculation (from field) ... ', end='')
srwl.CalcPowDenSR(stkP, eBeam, 0, magFldCnt, arPrecP)
print('done')

#**********************Saving results


# open output file
fs = open(outFile, 'wb')
header="#F "+outFile+" \n"
fs.write(header)

# dump inputs
header=''
header = header+ '#U ElectronEnergy = ' + repr(ElectronEnergy) + '\n'
header = header+ '#U ElectronCurrent = ' + repr(ElectronCurrent ) + '\n'
header = header+ '#U ElectronBeamSizeH = ' + repr(ElectronBeamSizeH ) + '\n'
header = header+ '#U ElectronBeamSizeV = ' + repr(ElectronBeamSizeV ) + '\n'
header = header+ '#U ElectronBeamDivergenceH = ' + repr(ElectronBeamDivergenceH ) + '\n'
header = header+ '#U ElectronBeamDivergenceV = ' + repr(ElectronBeamDivergenceV ) + '\n'
header = header+ '#U PeriodID = ' + repr(PeriodID ) + '\n'
header = header+ '#U N = ' + repr(N ) + '\n'
header = header+ '#U Kv = ' + repr(Kv ) + '\n'
header = header+ '#U PhotonEnergyMin = ' + repr(PhotonEnergyMin ) + '\n'
header = header+ '#U PhotonEnergyMax = ' + repr(PhotonEnergyMax ) + '\n'
header = header+ '#U PhotonEnergyPoints = ' + repr(PhotonEnergyPoints ) + '\n'
header = header+ '#U d = ' + repr(d ) + '\n'
header = header+ '#U gapH = ' + repr(gapH ) + '\n'
header = header+ '#U gapV = ' + repr(gapV ) + '\n'
header = header+ '#U B0 = ' + repr(B0 ) + '\n'
fs.write(header)


#
# write trajectory
#
header="\n#S 1  electron trajectory\n"
fs.write(header)
header="#N 7 \n#L ct[m]  X[m]  BetaX[rad]  Y[m]  BetaY[rad]  Z[m]  BetaZ[m]\n"
fs.write(header)
ctStep = 0
if partTraj.np > 0:
    ctStep = (partTraj.ctEnd - partTraj.ctStart)/(partTraj.np - 1)
ct = partTraj.ctStart
for i in range(partTraj.np):
    fs.write(str(ct) + '  ' + repr(partTraj.arX[i]) + '  ' + repr(partTraj.arXp[i]) + '  ' + repr(partTraj.arY[i]) + '  ' + repr(partTraj.arYp[i]) + '  ' + repr(partTraj.arZ[i]) + '  ' + repr(partTraj.arZp[i]) + '\n') 
    ct += ctStep

#
# write power density (mesh scan)
#
header="\n#S 2  power density W/mm2\n"
fs.write(header)
header="#N 3 \n#L V[mm]  H[mm]  PowerDensity[W/mm^2] \n"
fs.write(header)

for i in range(stkP.mesh.nx):  
    for j in range(stkP.mesh.ny): 
        xx = stkP.mesh.xStart + i*(stkP.mesh.xFin-stkP.mesh.xStart)/(stkP.mesh.nx-1)
        yy = stkP.mesh.yStart + j*(stkP.mesh.yFin-stkP.mesh.yStart)/(stkP.mesh.ny-1)
        ij = i*stkP.mesh.nx + j
        fs.write(repr(xx) + ' ' + repr(yy) + ' ' + repr(stkP.arS[ij]) + '\n')

#
# write flux
#
codata_ec = 1.602176565e-19 # electron charge in Coulomb http://physics.nist.gov/cgi-bin/cuu/Value?e
header="\n#S 3  flux\n"
fs.write(header)
header="#N 5 \n#L PhotonEnergy[eV]  Flux[phot/s/0.1%bw]  Flux[Phot/s/eV]  Power[Watts/0.1%bw]  Power[Watts/eV] \n"
fs.write(header)
for i in range(stkF.mesh.ne): 
    ener = stkF.mesh.eStart+i*(stkF.mesh.eFin-stkF.mesh.eStart)/(stkF.mesh.ne-1)
    flx1 = stkF.arS[i]       #  phot/sec/0.1%bw
    flx2 = flx1/ener*1e3     #  phot/sec/eV
    pwr1 = flx1*codata_ec/(0.1e-2)  #  W/0.1%bw
    pwr2 = flx1*ener*codata_ec      #  W/eV
    fs.write(' ' + repr(ener) + '   ' + 
        repr(flx1) + '   ' + repr(flx2) + '   ' + repr(pwr1) + '   ' + repr(pwr2) +
        '\n')
fs.close()

print ("File written to disk: ",outFile)
print ("Done")

