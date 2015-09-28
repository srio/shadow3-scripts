#
# script for running SRW to create a SHADOW source
#
import json
import numpy
import srwlib as sl
import Shadow as sd
import array
import cPickle as pickle
import sys


def ElectronBeam(x=0., y=0., z=0., xp=0., yp=0., e=6.04, Iavg=0.2, sigX=345e-6*1.e-20, sigY=23e-6*1.e-20, mixX=0.0, mixY=0.0, sigXp=4.e-9*1.e-20/345e-6, sigYp=4.e-11*1.e-20/23e-6, sigE = 1.e-4):
  el_rest = 0.51099890221e-03
  eBeam = sl.SRWLPartBeam()
  eBeam.Iavg = Iavg
  eBeam.partStatMom1.x     =  x
  eBeam.partStatMom1.y     =  y
  eBeam.partStatMom1.z     =  z
  eBeam.partStatMom1.xp    =  xp
  eBeam.partStatMom1.yp    =  yp
  eBeam.partStatMom1.gamma =  e/el_rest
  eBeam.partStatMom1.relE0 =  1.0
  eBeam.partStatMom1.nq    = -1
  eBeam.arStatMom2[ 0] = sigX**2  #from here it is not necessary for Single Electron calculation, obviously....
  eBeam.arStatMom2[ 1] = mixX
  eBeam.arStatMom2[ 2] = sigXp**2
  eBeam.arStatMom2[ 3] = sigY**2
  eBeam.arStatMom2[ 4] = mixY
  eBeam.arStatMom2[ 5] = sigYp**2
  eBeam.arStatMom2[10] = sigE**2
  return eBeam


def DriftElectronBeam(eBeam, und ):
  if isinstance(und, float):
    length = und
  elif isinstance(und, sl.SRWLMagFldU):    # Always defined in (0., 0., 0.) move the electron beam before the magnetic field.
    length = 0.0-0.55*und.nPer*und.per-eBeam.partStatMom1.z
  elif isinstance(und, sl.SRWLMagFldC): 
    if isinstance(und.arMagFld[0], sl.SRWLMagFldU):
      length = und.arZc[0]-0.55*und.arMagFld[0].nPer*und.arMagFld[0].per-eBeam.partStatMom1.z
    else: raise NameError
  else: raise NameError
  eBeam.partStatMom1.z += length
  eBeam.arStatMom2[0]  += 2*length*eBeam.arStatMom2[1]+length**2*eBeam.arStatMom2[2]
  eBeam.arStatMom2[1]  += length*eBeam.arStatMom2[2]
  eBeam.arStatMom2[3]  += 2*length*eBeam.arStatMom2[4]+length**2*eBeam.arStatMom2[5]
  eBeam.arStatMom2[4]  += length*eBeam.arStatMom2[5]
  eBeam.moved = length
  return eBeam  


def SimpleUndulator(nPer=72, per=0.0228, B=0.120215, n=1, h_or_v='v'):
  harmB = sl.SRWLMagFldH(n, h_or_v, B)
  und = sl.SRWLMagFldU([harmB], per, nPer)
  return und


def Undulator(nPer=72, per=0.0228, B=[0.120215], n=[1], h_or_v=['v']):
  assert (len(B)==len(n)), "Wrong length of input arrays"
  assert (len(B)==len(h_or_v)), "Wrong length of input arrays"
  harms = [ sl.SRWLMagFldH(n[i], h_or_v[i], B[i]) for i in range(len(B)) ]
  und = sl.SRWLMagFldU(harms, per, nPer)
  return und

  
def Undulators(und, Xc, Yc, Zc):#for the moment only one works
  cnt = sl.SRWLMagFldC([und], array.array('d', [Xc]), array.array('d', [Yc]), array.array('d', [Zc]))
  return cnt


#def SrwSESource(eBeam, cnt, mesh=sl.SRWLRadMesh(12000., 16000., 101, -15.e-6*50*3, 15e-6*50*3, 61, -15e-6*50*3, 15e-6*50*3, 61, 50.),  params=[1, 0.01, 0., 0., 20000, 1, 0]):
def SrwSESource(eBeam, cnt, mesh=sl.SRWLRadMesh(14718.4-1, 14718.4+1., 101, -15.e-6*50*3, 15e-6*50*3, 61, -15e-6*50*3, 15e-6*50*3, 61, 50.),  params=[1, 0.01, 0., 0., 20000, 1, 0]):
  wfr = sl.SRWLWfr()
  wfr.mesh = mesh
  wfr.partBeam = eBeam
  wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
  eBeam = DriftElectronBeam(eBeam, cnt)
  sl.srwl.CalcElecFieldSR(wfr, 0, cnt, params)
  stk = sl.SRWLStokes()
  stk.mesh = mesh
  stk.allocate(mesh.ne, mesh.nx, mesh.ny)
  eBeam = DriftElectronBeam(eBeam, -eBeam.moved)
  wfr.calc_stokes(stk)
  return stk, eBeam


def SrwMESource(eBeam, und, mesh=sl.SRWLRadMesh(14718.4, 14718.4, 1, -15.e-6*50, 15e-6*50, 81, -15e-6*50, 15e-6*50, 81, 50.),  params=[1, 9, 1.5, 1.5, 2]):
#def SrwMESource(eBeam, und, mesh=sl.SRWLRadMesh(1000., 21000., 10001, -15.e-6*50, 15e-6*50, 1, -15e-6*50, 15e-6*50, 1, 50.),  params=[1, 21, 1.5, 1.5, 1]):
  stk = sl.SRWLStokes()
  stk.mesh = mesh
  stk.allocate(mesh.ne, mesh.nx, mesh.ny)
  sl.srwl.CalcStokesUR(stk, eBeam, und, params)
  return stk, eBeam


def save(stk, eBeam, fname="SrwStokes"):
  pickle.dump( stk, open(fname+"_stk.dat", "wb") )
  pickle.dump( eBeam, open(fname+"_ebeam.dat", "wb") )


def Stokes0ToSpec(stk, fname="srw_xshundul.spec"):
  Shape = (4,stk.mesh.ny,stk.mesh.nx,stk.mesh.ne)
  data = numpy.ndarray(buffer=stk.arS, shape=Shape,dtype=stk.arS.typecode)
  data0 = data[0]
  x = numpy.linspace(stk.mesh.xStart,stk.mesh.xFin,stk.mesh.nx)
  y = numpy.linspace(stk.mesh.yStart,stk.mesh.yFin,stk.mesh.ny)
  e = numpy.linspace(stk.mesh.eStart,stk.mesh.eFin,stk.mesh.ne)
  f = open(fname,"w")
  for k in range(len(e)):
    f.write("#S %d intensity E= %f\n"%(k+1,e[k]))
    f.write("#N 3\n")
    f.write("#L X[m]  Y[m]  Intensity\n")
    for i in range(len(x)):
      for j in range(len(y)):
        f.write( "%e   %e   %e\n"%(x[i], y[j], data0[j,i,k]))
  f.close()
  sys.stdout.write('  file written: srw_xshundul.spec\n')
  
if __name__=="__main__":


  
  #
  # read inputs from a file created by ShadowVUI ----------------------------
  #
  inFileTxt = "xshundul.json"
  with open(inFileTxt, mode='r') as f1:
      h = json.load(f1)
  
  # list all non-empty keywords
  print "-----------------------------------------------------"
  for i,j in h.items():
      if (j != None):
          print "%s = %s" % (i,j)
  print "-----------------------------------------------------"
  print "k: ",h['K']
  
  #lambdau = 0.022888
  #k = 0.256
  #e_energy = 6.04
  #nperiods = 72
  #emin = 14705.0
  #emax = 14711.0
  #intensity = 0.2
  #maxangle = 0.015
  #sx = 345e-4
  #sz = 23e-4
  #ex = 4e-7
  #ez = 4e-9
  #nrays = 5000
  
  lambdau = h['LAMBDAU']
  k = h['K']
  e_energy = h['E_ENERGY']
  nperiods = h['NPERIODS']
  emin = h['EMIN']
  emax = h['EMAX']
  intensity = h['INTENSITY']
  maxangle = h['MAXANGLE']
  sx = h['SX']
  sz = h['SZ']
  ex = h['EX']
  ez = h['EZ']
  nrays = h['NRAYS']
  
  
  print "lambdau = ",lambdau
  print "k = ",k
  print "e_energy = ",e_energy
  print "nperiods = ",nperiods
  print "emin = ",emin
  print "emax = ",emax
  print "intensity = ",intensity
  print "maxangle = ",maxangle
  print "sx = ",sx
  print "sz = ",sz
  print "ex = ",ex
  print "ez = ",ez
  print "nrays = ",nrays
  print "emin = ",emin
  print "emax = ",emax
  
  
  #
  # define additional parameters needed by SRW
  #
  B = k/93.4/lambdau
  slit_distance = 50.0
  method = "SE" # single-electron  "ME" multi-electron
  nx = 51
  nz = 51
  sE = 1e-9 # 0.89e-3

  #
  # prepare inputs
  #

  # convert cm to m
  sx *= 1.0e-2 
  sz *= 1.0e-2
  ex *= 1.0e-2
  ez *= 1.0e-2
  
  sxp = ex/sx
  szp = ez/sz
  xxp = 0.0
  zzp = 0.0
  

  estep = 1.0
  # get a reasonable number of rays (estep = 1eV or 10 eV or 100 eV etc.)
  ne = 1+(emax-emin)/estep
  while ne > 100:
      ne = ne/10
  ne = int(ne)
  print "ne = ",ne
  #paramSE = [1, 0.01, 0, 0, 20000, 1, 0]
  paramSE = [1, 0.01, 0, 0, 50000, 1, 0]
  paramME = [1, 9, 1.5, 1.5, 2]

  # 
  # 
  if nx==1 and nz==1: paramME[4] = 1
  params = paramSE if method=="SE" else paramME

  
  slit_xmin = -maxangle*1.0e-3*slit_distance  
  slit_xmax =  maxangle*1.0e-3*slit_distance
  slit_zmin = -maxangle*1.0e-3*slit_distance
  slit_zmax =  maxangle*1.0e-3*slit_distance

  #
  # calculations
  #
  print("nperiods: %d, lambdau: %f, B: %f)"%(nperiods,lambdau,B))
  und = SimpleUndulator(nperiods,lambdau,B)
  print("e=%f,Iavg=%f,sigX=%f,sigY=%f,mixX=%f,mixY=%f,sigXp=%f,sigYp=%f,sigE=%f"%(e_energy,intensity,sx,sz,xxp,zzp,sxp,szp,sE) )
  eBeam = ElectronBeam(e=e_energy,Iavg=intensity,sigX=sx,sigY=sz,mixX=xxp,mixY=zzp,sigXp=sxp,sigYp=szp,sigE=sE)
  cnt = Undulators(und, 0., 0., 0.)
  sys.stdout.write('  calculating SE...'); sys.stdout.flush()
  print("emin=%f,emax=%f,ne=%d,slit_xmin=%f,slit_xmax=%f,nx=%d,slit_zmin=%f,slit_zmax=%f,nz=%d,slit_distance=%f"%(emin,emax,ne,slit_xmin,slit_xmax,nx,slit_zmin,slit_zmax,nz,slit_distance) )
  mesh = sl.SRWLRadMesh(emin,emax,ne,slit_xmin,slit_xmax,nx,slit_zmin,slit_zmax,nz,slit_distance)
  if (method == 'SE'):
      print ("Calculating SE...")
      stkSE, eBeam = SrwSESource(eBeam, cnt, mesh, params)
      sys.stdout.write('  done\n')
      sys.stdout.write('  saving SE Stokes...'); sys.stdout.flush()
      Stokes0ToSpec(stkSE)
      stk = stkSE
  else:
      print ("Calculating ME...")
      stkME, eBeam = SrwMESource(eBeam, und) # cnt, mesh, params)
      sys.stdout.write('  done\n')
      sys.stdout.write('  saving SE Stokes...'); sys.stdout.flush()
      Stokes0ToSpec(stkME)
      stk = stkME
  
  beam, param = sd.ShadowSrw.genShadowBeam( (stk,eBeam,), N=nrays, method=method, energy=None, lim=None, canted=None, distance=slit_distance)
  beam.write("begin.dat")
  sys.stdout.write('  file written: begin.dat\n')
  
  #eBeam = ElectronBeam()
  #print '  calculating ME'
  #stkME8, eBeam = SrwMESource(eBeam, und)
  #print '  done'
  #save(stkME, eBeam, "SrwStokesME")


