import numpy as np
from silx import config
from silx import sx
from silx.gui.plot.ScatterView import ScatterView
from PyQt5.QtWidgets import QApplication
import sys


if __name__ == "__main__":
    app = QApplication(sys.argv)

    samps = 1000
    levels = 100
    fwhm = 50
    cen = 500

    x = np.random.rand(samps)*levels
    y = np.random.rand(samps)*levels
    v = np.random.rand(samps)*levels
    x0, y0 = x[cen], y[cen]
    #2D Gaussian kernel to apply transparency
    a = np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

    config.DEFAULT_COLORMAP_NAME = 'viridis'

    #without transparency
    s1 = ScatterView()
    s1.setData(x, y, v)
    s1.show()
    #with transparency
    s2 = ScatterView()
    s2.setData(x, y, v, alpha=a)
    s2.show()

    app.exec()