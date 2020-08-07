import numpy
from srxraylib.plot.gol import plot, set_qt

set_qt()


y = numpy.linspace(-0.5,0.5,500)
period = 0.1
x = numpy.sin(2 * numpy.pi * y / 0.1)


filename = "/Users/srio/Oasys/Bsin.txt"
f = open(filename, 'w')
for i in range(x.size):
    f.write("%g %g\n" % (y[i], x[i]))
f.close()
print("File %s written to disk." % filename)

plot(y,x)