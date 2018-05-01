__author__ = 'srio'


import matplotlib.pyplot as plt
import numpy as np

from pynx.wavefront import *
from pynx.wavefront.cpu_operator import *  # Force using CPU operator
from pynx.utils.plot_utils import cm_phase


from srxraylib.plot.gol import plot_image
from numpy.testing import assert_almost_equal



pixel_size=10e-9
w = Wavefront(d='ascent', pixel_size=pixel_size, wavelength=1.5e-10)
w = ImshowAbs() * CircularMask(2e-6) * w

# w1 = ImshowAbs(fig_num=2) * MagnifyNearField(1.7, verbose=True) * w.copy()
# print("w1.z",w1.z)
# w2 = ImshowAbs(fig_num=3) * PropagateNearField(dz=w1.z, verbose=True) * w.copy()
# x0,x1 = plt.xlim()
# plt.figure(2)
# plt.xlim(x0,x1)
# plt.ylim(x0,x1)


w1 = ImshowAngle(fig_num=4) * MagnifyNearField(1.7, verbose=True) * w.copy()
x,y = w1.get_x_y()
z = w1.get(shift=True)
print(x.shape,y.shape,z.shape)
w2 = ImshowAngle(fig_num=5) * PropagateNearField(dz=w1.z, verbose=True) * w.copy()
x0,x1 = plt.xlim()
plt.figure(4)
plt.xlim(x0,x1)
plt.ylim(x0,x1)


x,y = w.get_x_y()
print("xmin: %f,xmax: %f,ymin: %f,ymax: %f"%(x.min(),x.max(),y.min(),y.max()))
X = np.linspace(y.min(),y.max(),y.size)
Y = np.linspace(y.min(),x.max(),x.size)
z = (w.get(shift=True))[0,:,:]
Z = np.flip(np.swapaxes(z,0,1),1)
# plot_image(np.abs(z1),1e6*X,1e6*Y,title="ORIG",show=False)



#
# wofry propagation
#

from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
from wofry.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
from wofry.beamline.optical_elements.ideal_elements.screen import WOScreen
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.beamline_element import BeamlineElement



ww = GenericWavefront2D.initialize_wavefront_from_arrays(X,Y,Z,wavelength=1.5e-10)
propagator = PropagationManager.Instance()
propagator.add_propagator(FresnelZoomXY2D())


screen = WOScreen(name="PIRRONE")
coordinates = ElementCoordinates(p = 0.0, q=w1.z)

propagation_elements = PropagationElements()
propagation_elements.add_beamline_element(BeamlineElement(optical_element=screen,
                                                          coordinates=coordinates))

parameters = PropagationParameters(wavefront=ww,
                                   propagation_elements=propagation_elements)
parameters.set_additional_parameters("shift_half_pixel", 1)
parameters.set_additional_parameters("magnification_x", 1.7)
parameters.set_additional_parameters("magnification_y", 1.7)


output_ww = propagator.do_propagation(propagation_parameters=parameters,
                                        handler_name=FresnelZoomXY2D.HANDLER_NAME)
plot_image(np.angle(output_ww.get_complex_amplitude()),
           1e6*output_ww.get_coordinate_x(),
           1e6*output_ww.get_coordinate_y(),
           title="WOFRY PHASE",cmap=cm_phase,show=0)

plt.xlim(x0,x1)
plt.ylim(x0,x1)



z1 = (w1.get(shift=True))[0,:,:]
Z1 = np.flip(np.swapaxes(z1,0,1),1)
plot_image(np.angle(Z1),
           1e6*output_ww.get_coordinate_x(),
           1e6*output_ww.get_coordinate_y(),
           title="PYNX PHASE",cmap=cm_phase,show=0)

plt.xlim(x0,x1)
plt.ylim(x0,x1)

plt.show()
