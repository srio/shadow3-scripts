

# 6,CPMU18a,UpStream,18.0,2000.00,  6.000, 2.6600, 0.0000, 0.8720, 1.0826, 1.0000, 1.1456, ,
# A = [2.6600, 0.0, 0.8720, 1.0826, 1.0000, 1.1456]
# id_period_mm = 18.0


# 19,W150b,Middle,150.0,1030.00, 26.500, 2.6479, 0.0000, 0.8699, 1.2400, 1.0000, 1.1300, ,
A = [2.6479, 0.0000, 0.8699, 1.2400, 1.0000, 1.1300]
id_period_mm = 150.0
min_gap = 26.500
# self.KY = (93.4 * 0.01 * self.PERIOD * (
            # 2.3333 * numpy.exp(-0.02473 * self.GAP) + 1.189 * numpy.exp(-0.059691 * self.GAP)))

# 1,U35a,UpStream,35.0,1600.00, 11.000, 1.9562, 1.0131, , , , , ,
# A = [1.9562, 1.0131, 0, 0, 0, 0]
#
# id_period_mm = 16.0


import numpy
import scipy.constants as codata

# A = [3.29, 1.095, 0, 0, 0, 0] # from Reine+Joel report CPMU value


gap_mm = numpy.linspace(min_gap, 2*min_gap, 100)


Bmax = numpy.zeros_like(gap_mm)

Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * gap_mm / id_period_mm)
Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * gap_mm / id_period_mm)
Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * gap_mm / id_period_mm)

# Bmax += A[0] * numpy.exp(-numpy.pi * A[1] * gap_mm / id_period_mm)
# Bmax += A[2] * numpy.exp(-numpy.pi * A[3] * gap_mm / id_period_mm)
# Bmax += A[4] * numpy.exp(-numpy.pi * A[5] * gap_mm / id_period_mm)


K = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)
# self.KY = (93.4 * 0.01 * self.PERIOD * (
            # 2.3333 * numpy.exp(-0.02473 * self.GAP) + 1.189 * numpy.exp(-0.059691 * self.GAP)))
K2 = (93.4 * 0.01 * (id_period_mm/10) * (
            2.3333 * numpy.exp(-0.02473 * gap_mm) + 1.189 * numpy.exp(-0.059691 * gap_mm)))

print(Bmax, K)

from srxraylib.plot.gol import plot
plot(gap_mm, K, gap_mm, K2, xtitle="gap [mm]", ytitle="K", legend=['now','gwen'])

