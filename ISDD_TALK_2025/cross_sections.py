
#
# script to make the calculations (created by XOPPY:crosssec)
#
from xoppylib.scattering_functions.xoppy_calc_crosssec import xoppy_calc_crosssec
import xraylib
from dabax.dabax_xraylib import DabaxXraylib

from srxraylib.plot.gol import plot


"""
    :param calculate:
            0: total cross section
            1: photoelectric cross section
            2: rayleigh cross serction
            3: compton cross section
            4: total minus raileigh cross section
"""


#
# Diamond
#
title="Diamond (Z=6)"
OUT_DICT = []
for calculate in [2,1,3,0]:
    out_dict =  xoppy_calc_crosssec(
        descriptor   = "C",
        density      = 3.5,
        MAT_FLAG     = 0,
        CALCULATE    = calculate,
        GRID         = 1,
        GRIDSTART    = 200.0,
        GRIDEND      = 500000.0,
        GRIDN        = 2000,
        UNIT         = 0,
        DUMP_TO_FILE = 0,
        FILE_NAME    = "CrossSec.dat",
        material_constants_library = xraylib,
        )
    OUT_DICT.append(out_dict)


plot(
    OUT_DICT[0]["data"][0,:],OUT_DICT[0]["data"][-1,:],
    OUT_DICT[1]["data"][0,:],OUT_DICT[1]["data"][-1,:],
    OUT_DICT[2]["data"][0,:],OUT_DICT[2]["data"][-1,:],
    OUT_DICT[3]["data"][0,:],OUT_DICT[3]["data"][-1,:],
    xtitle=OUT_DICT[0]["labels"][0],
    ytitle="",
    legend=[OUT_DICT[0]["labels"][1],OUT_DICT[1]["labels"][1],OUT_DICT[2]["labels"][1],OUT_DICT[3]["labels"][1]],
    title=title,
    linestyle=[None,None,None,'--'], figsize=(8,6),
    xlog=True,ylog=True,show=0)


#
# Copper
#
title="Copper (Z=29)"
OUT_DICT = []
for calculate in [2,1,3,0]:
    out_dict =  xoppy_calc_crosssec(
        descriptor   = "Cu",
        density      = 8.96,
        MAT_FLAG     = 0,
        CALCULATE    = calculate,
        GRID         = 1,
        GRIDSTART    = 200.0,
        GRIDEND      = 500000.0,
        GRIDN        = 2000,
        UNIT         = 0,
        DUMP_TO_FILE = 0,
        FILE_NAME    = "CrossSec.dat",
        material_constants_library = xraylib,
        )
    OUT_DICT.append(out_dict)


plot(
    OUT_DICT[0]["data"][0,:],OUT_DICT[0]["data"][-1,:],
    OUT_DICT[1]["data"][0,:],OUT_DICT[1]["data"][-1,:],
    OUT_DICT[2]["data"][0,:],OUT_DICT[2]["data"][-1,:],
    OUT_DICT[3]["data"][0,:],OUT_DICT[3]["data"][-1,:],
    xtitle=OUT_DICT[0]["labels"][0],
    ytitle="",
    legend=[OUT_DICT[0]["labels"][1],OUT_DICT[1]["labels"][1],OUT_DICT[2]["labels"][1],OUT_DICT[3]["labels"][1]],
    title=title,
    linestyle=[None,None,None,'--'], figsize=(8,6),
    xlog=True,ylog=True,show=True)





