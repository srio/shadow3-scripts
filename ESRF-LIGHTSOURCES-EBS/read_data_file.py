import numpy


filename = 'IDstatus_2020-10-23.dat'

a = numpy.loadtxt(filename, skiprows=1, dtype=str)

print(a.shape)

straight_section = a[:,0].astype(int)
id_name = a[:,1]
id_period = 1e-3 * a[:,2].astype(float)
id_length = 1e-3 * a[:,3].astype(float)
id_minimum_gap_mm = a[:,4].astype(float)
a1 = a[:,5].astype(float)
a2 = a[:,6].astype(float)
a3 = a[:,7].astype(float)
a4 = a[:,8].astype(float)
a5 = a[:,9].astype(float)
a6 = a[:,10].astype(float)
a7 = a[:,11].astype(float)

print(a1)
a0 = id_minimum_gap_mm
Bmax = numpy.zeros_like(a1)
Bmax += a1 * numpy.exp( -1 * numpy.pi * (id_minimum_gap_mm - a0 ))
Bmax += a2 * numpy.exp( -2 * numpy.pi * (id_minimum_gap_mm - a0 ))
Bmax += a3 * numpy.exp( -3 * numpy.pi * (id_minimum_gap_mm - a0 ))
Bmax += a4 * numpy.exp( -4 * numpy.pi * (id_minimum_gap_mm - a0 ))
Bmax += a5 * numpy.exp( -5 * numpy.pi * (id_minimum_gap_mm - a0 ))
Bmax += a6 * numpy.exp( -6 * numpy.pi * (id_minimum_gap_mm - a0 ))
Bmax += a7 * numpy.exp( -7 * numpy.pi * (id_minimum_gap_mm - a0 ))

print(Bmax)



