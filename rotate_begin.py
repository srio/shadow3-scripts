import Shadow
import numpy

# shift and rotate source

# rotate ray directions an angle theta around the Z axis (vertical)


theta = 0.0062
y0 = -55.0
x0 = 0.5902 ; 0.579

print("Rotation angle is: %f degrees"%(theta*180/numpy.pi))
a = Shadow.Beam()
a.load("/users/srio/Working/rt/ESRF-new-lattice/begin.dat")
print(a.rays.shape)

x =  a.rays[:,0]  # column 1 x
y =  a.rays[:,1]  # column 2 y
xp = a.rays[:,3]  # column 4 xprime
yp = a.rays[:,4]  # column 5 yprime

# help,x,yxp,yp
# 
# 
xp_new =  xp*numpy.cos(theta) + yp*numpy.sin(theta)
yp_new = -xp*numpy.sin(theta) + yp*numpy.cos(theta)


y = y + y0

x_new =  x*numpy.cos(theta) + y*numpy.sin(theta)
y_new = -x*numpy.sin(theta) + y*numpy.cos(theta)

y_new = y_new - y0


a.rays[:,0] = x_new
a.rays[:,1] = y_new
a.rays[:,3] = xp_new
a.rays[:,4] = yp_new
a.rays[:,0] = a.rays[:,0] + x0

a.write("begin2.dat")

Shadow.ShadowTools.plotxy(a,2,1,nbins=201)
#
# xsh_plotxy,'begin.dat',2,1,NoLost=0,Histo=2,NBins=151,CalFWHM=1,CCol=2,title='theta='+strcompress(theta)
# xsh_plotxy,'star.01',1,3,Retrace=0.0000,NoLost=0,Histo=2,NBins=151,CalFWHM=1,CCol=2,CFile='begin.dat'
# xsh_plotxy,'mirr.01',2,1,NoLost=0,Histo=2,NBins=151,CalFWHM=1,CCol=2,CFile='begin.dat'
# 
