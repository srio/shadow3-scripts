import numpy as np
import barc4plots as b4pt
from barc4plots.barc4plots import Image2Plot
from barc4plots.barc4plots import plot_2D_cuts, plot_2D, ESRF_colors_2D
from srxraylib.plot.gol import set_qt, plot_image_with_histograms, plot_show

set_qt()

nx = 100
ny = 100
image = np.zeros((ny, nx))
x = np.linspace(-10, 10, nx)
y = np.linspace(-1, 1, ny)
X = np.outer(x, np.ones_like(y))
Y = np.outer(np.ones_like(x), y)
image = X * Y

# plot_container = Image2Plot(image, x, y)
# # plot_container.ax_limits = [-4., 3., -1.1, 1.1]
# # plot_container.legends = ['2D plot testing - intensity cuts', 'x-axis (au)', 'y-axis (au)']
# plot_container.AspectRatio = True
# plot_container.ColorScheme = 8  # see ESRF_colours() and ESRF_colours_2D()
# # plot_container.plt_limits = [0, 1]
# # plot_container.Scale = 0  # 2D plot: 0 - linear; 1 - log10; 2: Gamma = 0.25.
# # plot_container.sort_class()
#
# plot_2D(plot_container, Crop=False, ROI='r', Scale=False, Enable=False, Silent=False, isphase=False,
#         m=6.4, n=4.8, dpi=300)
#
# plot_2D_cuts(plot_container, Enable=True, Silent=False, isphase=False, m=7.766563146, n=4.8, dpi=300, x=None, y=None)


# new

image = Image2Plot(image, x, y)
image.legends = ['ideal case', '($\mu$m)', '(mm)']
image.AspectRatio = False
# # image.ax_limits = [-25, 25, -2.5, 2.5]
image.Scale = 0
image.sort_class()
plot_2D_cuts(image, Enable=True, Silent=False, m=7.766563146, n=4.8)



# if plot_container.AspectRatio:
#     aspect_ratio = 'equal'
# else:
#     aspect_ratio = 'auto'
#
# plot_image_with_histograms(image, x, y, cmap=ESRF_colors_2D(8), aspect_ratio=aspect_ratio, figsize=(10,1.5), show=0,
#                            use_profiles_instead_histograms=True)
#
#
# plot_show()