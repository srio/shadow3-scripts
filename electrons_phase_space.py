#
# random sampling of the electron phase space (tests...)
#

import numpy 
import matplotlib.pyplot as plt

# H
meanH = [0,0]
covH = [[10,50],[50,100]] # diagonal covariance, points lie on x or y-axis

x,xp = numpy.random.multivariate_normal(meanH,covH,5000).T


plt.plot(x,xp,'.')
plt.axis('equal')
plt.show()
