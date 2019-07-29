import numpy
from srxraylib.metrology.dabam import dabam
from srxraylib.plot.gol import plot


if __name__ == "__main__":
    filename = "C:\\Users\\Manuel\\OASYS1.2\\shadow3-scripts\\DABAM_DANIELE\\TXI_Flat1.dat"

    a  = numpy.loadtxt(filename,skiprows=1)

    print(a.shape)

    plot(a[:,0],a[:,3])

    d = dabam()
    d.load()