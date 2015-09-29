import numpy as np
import matplotlib.pyplot as plt


# inputs (mean is zero) 
mean = [0,0]
sig1 = 1.0
sig2 = 2.0
rho = -0.75
Npoints = 5000


#covariance matrix
cov = np.array( [[sig1*sig1,rho*sig1*sig2],[rho*sig1*sig2,sig2*sig2]] )
print("\n\n input covariance matrix: ",cov)

#
# method 1: using np routine
#
x,y = np.random.multivariate_normal(mean,cov,Npoints).T
plt.plot(x,y,'x'); plt.axis('equal'); plt.title('method 1') ; plt.show()

#
# method 2: using Chalesky decomposition
#
from numpy import linalg as LA
U = LA.cholesky(cov)
XX = np.random.normal(size=2*Npoints)
XX.shape = (2,Npoints)
x,y = np.dot(U,XX) 
plt.plot(x,y,'x'); plt.axis('equal'); plt.title('method 2') ; plt.show()

#
# method 3: using eigenvvalue decomposition
#
#http://homepages.inf.ed.ac.uk/imurray2/code/matlab_octave_missing/mvnrnd.m

Lambda,E = LA.eig(cov)
LambdaM = [ [Lambda[0],0],[0,Lambda[1]]]
U = np.dot(np.sqrt(LambdaM),E.T)
U = U.T

XX = np.random.normal(size=2*Npoints)
XX.shape = (2,Npoints)
x,y = np.dot(U,XX) 
plt.plot(x,y,'x'); plt.axis('equal'); plt.title('method 3') ; plt.show()


# get statistics using np
a=np.cov( np.vstack((x,y))  )
print("\n\n reasult covariance matrix: ")
print(a)

#another way
X1 = np.array(x, ndmin=2, dtype=float)
Y1 = np.array(y, copy=False, ndmin=2, dtype=float)
axis = 0 
XX = np.concatenate((X1, Y1), axis)
N = X1.shape[1]
cc = (np.dot(XX, XX.T.conj()) / float(N-1) ).squeeze()
print(cc)
print("\n\n result covariance matrix (by hand): %",(cc))

# <x y>
mm = np.dot(x,y)/float(N)
print("\n\n cross term <x y> %f: "%(mm))
print("\n\n rho: %f "%(mm/np.sqrt(cc[0,0]*cc[1,1])))
print(mm)




