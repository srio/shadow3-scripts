

import numpy
from srxraylib.plot.gol import plot, plot_table, plot_show
from orangecontrib.xoppy.util.xoppy_xraylib_util import cross_calc_mix



energy = numpy.linspace(1.0, 300.0, 100) # in keV
density = 1.0
descriptor = "SiO2"
icalculate = 0
iunit = 3

descriptors = [
    "C",
    "SiO2",
    "B4C",
    "Al2O3",
    "PbO"
    ]

densities = [
    2.26,
    2.65,
    2.52,
    3.95,
    9.53,
    ]


calculate = [
            "0: total cross section",
            "1: photoelectric cross section",
            "2: rayleigh cross serction",
            "3: compton cross section",
            "4: total minus raileigh cross section",
            ]

unit = [
            "0: barn/atom (Cross Section calculation)",
            "1: cm^2 (Cross Section calculation)",
            "2: cm^2/g (Mass Attenuation Coefficient)",
            "3: cm^-1 (Linear Attenuation Coefficient)",
            ]
"""
    :param parse_or_nist: 0 for compound (default), 1 for name in the NIST compound list
    :param density: the material density in g/cm^3
"""




# #
# # cross sections
# #
# out = numpy.zeros((5,energy.size))
# for i,descriptor in enumerate(descriptors):
#     density = densities[i]
#     for icalc in range(5):
#         out[icalc,:] = cross_calc_mix(descriptor,1000*energy,calculate=icalc,unit=iunit,parse_or_nist=0,density=density)
#
#
#     plot_table(energy,out,xlog=True,ylog=True,
#                xtitle="Photon energy [kev]",ytitle=" %s"%(unit[iunit]),
#                title="Material: %s, density=%f"%(descriptor,density),
#                legend=calculate,legend_position=(0.98,0.98),
#                show=False,
#                )

# #
# # ratio compton/total
# #
# out = numpy.zeros((len(descriptors),energy.size))
# for i,descriptor in enumerate(descriptors):
#     density = densities[i]
#     out[i,:] = \
#         cross_calc_mix(descriptor,1000*energy,calculate=3,unit=iunit,parse_or_nist=0,density=density) \
#         / cross_calc_mix(descriptor,1000*energy,calculate=0,unit=iunit,parse_or_nist=0,density=density)
#
# plot_table(energy,out,xlog=True,ylog=True,
#            xtitle="Photon energy [kev]",ytitle=" %s"%(unit[iunit]),
#            title="Ratio Compton/Total",
#            legend=descriptors,legend_position=(0.98,0.98),
#            show=False,
#            )

#
# ratio compton/totalPbO
#
out = numpy.zeros((len(descriptors),energy.size))
for i,descriptor in enumerate(descriptors):
    density = densities[i]
    out[i,:] = \
        cross_calc_mix(descriptor,1000*energy,calculate=3,unit=iunit,parse_or_nist=0,density=density) \
        / \
        cross_calc_mix("SiO2",1000*energy,calculate=0,unit=iunit,parse_or_nist=0,density=2.65)

plot_table(energy,out,xlog=True,ylog=False,
           xtitle="Photon energy [kev]",ytitle=" %s"%(unit[iunit]),
           title="Ratio Compton/TotalSiO2",
           legend=descriptors,legend_position=(0.98,0.98),
           show=False,
           )


# #
# # ratio elastic/total
# #
# out = numpy.zeros((len(descriptors),energy.size))
# for i,descriptor in enumerate(descriptors):
#     density = densities[i]
#     out[i,:] = \
#         cross_calc_mix(descriptor,1000*energy,calculate=2,unit=iunit,parse_or_nist=0,density=density) \
#         / cross_calc_mix(descriptor,1000*energy,calculate=0,unit=iunit,parse_or_nist=0,density=density)
#
# plot_table(energy,out,xlog=True,ylog=True,
#            xtitle="Photon energy [kev]",ytitle=" %s"%(unit[iunit]),
#            title="Ratio Elastic/Total",
#            legend=descriptors,legend_position=(0.98,0.98),
#            show=False,
#            )

# #
# # ratio photoelectric/total
# #
# out = numpy.zeros((len(descriptors),energy.size))
# for i,descriptor in enumerate(descriptors):
#     density = densities[i]
#     out[i,:] = \
#         cross_calc_mix(descriptor,1000*energy,calculate=1,unit=iunit,parse_or_nist=0,density=density) \
#         / cross_calc_mix(descriptor,1000*energy,calculate=0,unit=iunit,parse_or_nist=0,density=density)
#
# plot_table(energy,out,xlog=True,ylog=True,
#            xtitle="Photon energy [kev]",ytitle=" %s"%(unit[iunit]),
#            title="Ratio photoelectric/Total",
#            legend=descriptors,legend_position=(0.98,0.98),
#            show=False,
#            )


plot_show()
