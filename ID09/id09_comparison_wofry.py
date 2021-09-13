import numpy
from srxraylib.plot.gol import plot


factor = 1.0 # 0.75

#
# 3 mrad
#
exp = numpy.loadtxt("/users/srio/Oasys/id09/measurements_2_998_mrad.csv",skiprows=1,delimiter=',')
sim = numpy.loadtxt("plot_energy_scan_3mrad.dat",skiprows=0)
hy = numpy.loadtxt("/users/srio/Oasys/id09/beam_2_998_15000.csv",skiprows=1,delimiter=',')

plot(
    exp[:,1], exp[:,2] / exp[:,2].mean(),
    factor * sim[:, 0], sim[:,1] / sim[:,1].mean() ,
    factor * sim[:, 0], sim[:,2] / sim[:,2].mean() / 5,
    hy[:,0], numpy.flip(hy[:,1] / hy[:,1].mean()),
    legend = ["experimental",
              "simulated Wofry1D (pink)",
              "simulated Wofry1D (monochromatic)",
              "Hybrid"],
    xrange=[-1200,1200], linestyle=[None,None,None,":"],
    color=["k","pink","r","b"], title="3.0 mrad", show=0,
    )


#
# 2.5 mrad
#
exp = numpy.loadtxt("/users/srio/Oasys/id09/measurements_2_498_mrad.csv",skiprows=1,delimiter=',')
sim = numpy.loadtxt("plot_energy_scan_2p5mrad.dat",skiprows=0)
hy = numpy.loadtxt("/users/srio/Oasys/id09/beam_2_498_15000.csv",skiprows=1,delimiter=',')

plot(
    exp[:,1], exp[:,2] / exp[:,2].mean(),
    factor * sim[:, 0], sim[:,1] / sim[:,1].mean() ,
    factor * sim[:, 0], sim[:,2] / sim[:,2].mean() / 5,
    hy[:,0], numpy.flip(hy[:,1] / hy[:,1].mean()),
    legend = ["experimental",
              "simulated Wofry1D (pink)",
              "simulated Wofry1D (monochromatic)",
              "Hybrid"],
    xrange=[-1200,1200],
    color=["k","pink","r","b"], linestyle=[None,None,None,":"], title="2.5 mrad",
    )

