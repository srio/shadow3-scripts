import numpy
from srxraylib.plot.gol import plot
srcalc = numpy.loadtxt("/users/srio/miniconda3/lib/python3.7/site-packages/xoppylib/data/reflect/Dia")
xrl = numpy.loadtxt("/scisoft/XRayOptics/OASYS1.2.1/miniconda3/lib/python3.7/site-packages/xoppylib/data/reflect/Dia")

plot(srcalc[:,0], srcalc[:,1],
     xrl[:,0], xrl[:,1],
srcalc[:,0], srcalc[:,2],
     xrl[:,0], xrl[:,2],
     xlog=1, ylog=1, legend=['Delta srcalc','Delta xraylib','Beta srcalc','Beta xraylib'], title='Be: delta')

# plot(srcalc[:,0], srcalc[:,2],
#      xrl[:,0], xrl[:,2],
#      xlog=1, ylog=1, legend=['srcalc','xraylib'], title='Be: beta')