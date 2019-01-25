# This example shows how to build colorbars without an attached mappable.
# https://matplotlib.org/examples/api/colorbar_only.html

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy

# Make a figure and axes with dimensions as desired.
fig = plt.figure(figsize=(8, 1))
# fig.tight_layout()
ax1 = fig.add_axes([0.05, 0.5, 0.9, 0.45])
# ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
# ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])
# ax3 = fig.add_axes([0.05, 0.15, 0.9, 0.15])

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = mpl.cm.hsv
norm = mpl.colors.Normalize(vmin=-numpy.pi, vmax=numpy.pi)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')

plt.xticks(fontsize=15)
ax1.xaxis.label.set_size(12)
cb1.set_label('phase [rad]')

plt.savefig("/tmp/colorbar.png")
print("File written to disk: /tmp/colorbar.png ")
plt.show()