import numpy
from srxraylib.plot.gol import plot
# sim = numpy.loadtxt("/users/srio/OASYS1.2/shadow3-scripts/ID09/results1/id09_3mrad_spectral_density_18489.dat",
#                     skiprows=4)
# current = numpy.flip(sim[:,1])
# plot(
#     sim[:, 0], current,
#     )


Energy = numpy.linspace(18000,22000,50)
Energy = numpy.linspace(18500,20500,100)

cumulated = 0.0
spectrum = numpy.zeros_like(Energy)
cf = numpy.zeros_like(Energy)


dir = "results_3mrad/"
dir = "results_2p5mrad/"

e0 = 20016.064
deltahalf = 1
for i,energy in enumerate(Energy):
    if energy > (e0-10000*deltahalf)  and energy < (e0+10000*deltahalf):
        if True: #try:
            sim = numpy.loadtxt("%s/id09_3mrad_spectral_density_%4d.dat" % (dir,energy),
                                skiprows=4)
            current = numpy.flip(sim[:,1])
            cumulated = cumulated + current
            spectrum[i] = current.sum() * (sim[1, 0] - sim[0, 0])


            cfi = numpy.loadtxt("%s/occupation_%4d.dat" % (dir, energy),
                                skiprows=4)
            cf[i] = cfi[0,1]
        # except:
        #     break


e00 = 20015
sim = numpy.loadtxt("%s/id09_3mrad_spectral_density_%4d.dat" % (dir, e00), skiprows=4)
current = numpy.flip(sim[:,1])

plot(
    sim[:, 0], current * 10,
    sim[:, 0], cumulated,
    xrange=[-1200,1200],
    legend=["monochromatic (E=%d eV) x 10" % e00, "pink"],
    color=["r","pink"]
    )

plot(Energy, spectrum)
plot(Energy, cf, title="CF")

filename = "plot_energy_scan_2p5mrad.dat"
f = open(filename, "w")
for i in range(cumulated.size):
    f.write("%g  %g  %g\n" % (sim[i, 0], cumulated[i], current[i]))
f.close()
print("File %s written to disk" % filename)