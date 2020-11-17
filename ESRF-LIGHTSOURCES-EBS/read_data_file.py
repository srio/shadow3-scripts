import numpy
import scipy.constants as codata

filename = 'IDstatus_2020-10-23.dat'

a = numpy.loadtxt(filename, skiprows=1, dtype=str)

print(a.shape)

straight_section = a[:,0].astype(int)
id_name = a[:,1]
id_period = 1e-3 * a[:,2].astype(float)
id_period_mm = a[:,2].astype(float)
id_length = 1e-3 * a[:,3].astype(float)
id_minimum_gap_mm = a[:,4].astype(float)
a0 = a[:,5].astype(float)
a1 = a[:,6].astype(float)
a2 = a[:,7].astype(float)
a3 = a[:,8].astype(float)
a4 = a[:,9].astype(float)
a5 = a[:,10].astype(float)
a6 = a[:,11].astype(float)


Bmax = numpy.zeros_like(a0)
Bmax += a1 * numpy.exp( -1 * numpy.pi * (id_minimum_gap_mm - a0 ) / id_period_mm)
Bmax += a2 * numpy.exp( -2 * numpy.pi * (id_minimum_gap_mm - a0 ) / id_period_mm)
Bmax += a3 * numpy.exp( -3 * numpy.pi * (id_minimum_gap_mm - a0 ) / id_period_mm)
Bmax += a4 * numpy.exp( -4 * numpy.pi * (id_minimum_gap_mm - a0 ) / id_period_mm)
Bmax += a5 * numpy.exp( -5 * numpy.pi * (id_minimum_gap_mm - a0 ) / id_period_mm)
Bmax += a6 * numpy.exp( -6 * numpy.pi * (id_minimum_gap_mm - a0 ) / id_period_mm)

Kmax = Bmax * id_period * codata.e /(2 * numpy.pi * codata.m_e * codata.c)

print("\n\n%5s  %10s  %15s %15s %15s" % ("sect", "name", "Gmin", "Bmax", "Kmax"))
print("%5s  %10s  %15s %15s %15s" % ("====", "====", "====", "====", "===="))
for i in range(Bmax.size):
    print("%5d  %10s  %15.3f %15.3f %15.3f" % (straight_section[i], id_name[i], id_minimum_gap_mm[i], Bmax[i], Kmax[i]))

out_dict = {}
out_dict["straight_section"] = straight_section.tolist()
out_dict["id_name"] = id_name.tolist()
out_dict["id_minimum_gap_mm"] = id_minimum_gap_mm.tolist()
out_dict["Bmax"] = Bmax.tolist()
out_dict["Kmax"] = Kmax.tolist()
out_dict["straight_section"] = straight_section.tolist()
out_dict["id_period"] = id_period.tolist()
out_dict["id_period_mm"] = id_period_mm.tolist()
out_dict["id_length"] = id_length.tolist()
out_dict["a0"] = a0.tolist()
out_dict["a1"] = a1.tolist()
out_dict["a2"] = a2.tolist()
out_dict["a3"] = a3.tolist()
out_dict["a4"] = a4.tolist()
out_dict["a5"] = a5.tolist()
out_dict["a6"] = a6.tolist()

import json



json_object = json.dumps(out_dict, indent = 4)
# print(json_object)

f = open('ebs_ids.json', 'w')
f.write(json_object)
f.close()
print("File written to disk: ebs_ids.json")


with open('ebs_ids.json') as json_file:
    data = json.load(json_file)

print(data)



