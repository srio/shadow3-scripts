#!/usr/bin/env python
# coding: utf-8

# In[78]:


import scipy.ndimage.filters as sc
import numpy as np
import matplotlib.pyplot as plt


# In[79]:


file='test.txt'
x,y,z=np.loadtxt(file, delimiter=';', unpack=True)


# In[80]:


xvals = np.unique(x)


# In[81]:


len(xvals)


# In[82]:


yvals=np.unique(y)
len(yvals)


# In[83]:


xvals


# In[84]:


yvals


# In[85]:


zvals=z.reshape(75,3178)


# In[86]:


zvals


# In[87]:


# get_ipython().run_line_magic('matplotlib', 'inline')


# In[88]:


imgplot = plt.imshow(zvals)


# In[89]:


fig, ax = plt.subplots()
ax.set_title('height profile')
ax.set_ylabel('height [m]')
ax.set_xlabel('position [m]')
ax.plot(xvals, zvals[37])  # Plot some data on the axes.


# In[93]:


#subtract best sphere/cylinder from line profile across cenre of mirror
quadfit=np.polyfit(xvals,zvals[37],2)
print(quadfit)
print('radius of curvature = %.2f m'% (1/(2*quadfit[0])))
height_err=zvals[37]-np.polyval(quadfit, xvals)
slope_err=(zvals[37,:-1]-zvals[37,1:])/(xvals[1]-xvals[0])
slope_err_xvals=(xvals[:-1]+xvals[1:])/2
#remove residual curvature
linfit=np.polyfit(slope_err_xvals,slope_err,1)
slope_err=slope_err-np.polyval(linfit, slope_err_xvals)
print(xvals[1]-xvals[0])
print('rms height err = %.2e m, pv height error = %.2e m'% (height_err.std(),height_err.max()-height_err.min()))
print('rms slope err = %.2e rad, pv slope error = %.2e rad'% (slope_err.std(),slope_err.max()-slope_err.min()))
#filter slope errors to reduce noise. Gaussian sigma 3 pixels equivalent psf of ~ 2mm fwhm for XFEL data.
filtered_slope_err=sc.gaussian_filter(slope_err,3) 
print('rms filtered slope err = %.2e rad, pv slope error = %.2e rad'% (filtered_slope_err.std(),filtered_slope_err.max()-filtered_slope_err.min()))


# In[52]:


fig, ax = plt.subplots()
ax.set_title('height error profile')
ax.set_ylabel('height [m]')
ax.set_xlabel('position [m]')
ax.plot(xvals, height_err)  # Plot some data on the axes.


# In[74]:


fig, ax = plt.subplots()
ax.set_title('slope error profile')
ax.set_ylabel('slope error [rad]')
ax.set_xlabel('position [m]')
ax.plot(slope_err_xvals, slope_err)  # Plot some data on the axes.


# In[56]:


print('rms height err = %.2e m, pv height error = %.2e m'% (height_err.std(),height_err.max()-height_err.min()))

#pv_height_err=height_err.max()-height_err.min()


# In[38]:


fig, ax = plt.subplots()
ax.set_title('height profile')
ax.set_ylabel('height [m]')
ax.set_xlabel('position [m]')
ax.plot(x, z)  # Plot some data on the axes.


# In[57]:


#fit whole surface with best cylinder
quadfit=np.polyfit(x,z,2)
print(quadfit)
print('radius of curvature = %.2f m'% (1/(2*quadfit[0])))
height_err_surf=z-np.polyval(quadfit, x)


# In[58]:


#2d array of surface height_error after best cylinder subtraction 
height_err_surf=height_err_surf.reshape(75,3178)


# In[59]:


imgplot = plt.imshow(height_err_surf)


# In[60]:


height_err_surf.std(), (height_err_surf.max()-height_err_surf.min())


# In[102]:


plt.pcolormesh(xvals,yvals,height_err_surf)




from srxraylib.plot.gol import plot_image
plot_image(height_err_surf, yvals, xvals, aspect='auto', show=0)

plot_image(z.reshape(75,3178), yvals, xvals, aspect='auto', show=0, title="RAW")


print(height_err_surf.shape, yvals.shape, xvals.shape)
# In[ ]:


plt.show()


# In[ ]:





# In[ ]:




