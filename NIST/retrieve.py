# Hubbell, J.H. and Seltzer, S.M. (2004), Tables of X-Ray Mass Attenuation Coefficients and Mass Energy-Absorption Coefficients (version 1.4). [Online] Available: http://physics.nist.gov/xaamdi [year, month day]. National Institute of Standards and Technology, Gaithersburg, MD.

import requests as req
import xraylib


txt = """#F MassEnergyAbsorption_NIST.dat
#C  This file belongs to the DABAX library. 
#UT Mass Energy Absorption Coefficients from NIST
#UD Mass Energy Absorption Coefficients from NIST
#UD 
#UD 
#UD Tables of X-Ray Mass Attenuation Coefficients and Mass Energy-Absorption
#UD Coefficients from 1 keV to 20 MeV for Elements Z = 1 to 92 and 
#UD 48 Additional Substances of Dosimetric Interest
#UD 
#UD Reference:
#UD Hubbell, J.H. and Seltzer, S.M. (2004), Tables of X-Ray Mass Attenuation Coefficients
#UD and Mass Energy-Absorption Coefficients (version 1.4). [Online] 
#UD Available: http://physics.nist.gov/xaamdi [2021-06-03].
#UD National Institute of Standards and Technology, Gaithersburg, MD.
"""

f = open("MassEnergyAbsorption_NIST.dat",'w')
f.write(txt)



for Z in range(1,93):
    print(">>> processing Z", Z)
    resp = req.get("https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z%02d.html" % Z)
    txt = resp.text
    # print(resp.text)

    f.write("\n#S %d %s\n#N 3\n#L  energy [eV]  mu/rho [cm2/g]  mu_energy/rho[cm2/g]\n" % (Z, xraylib.AtomicNumberToSymbol(Z)))

    txt_data = txt[txt.index("<PRE>"):txt.index("</PRE>")]
    # print(txt_data)

    txt_lines = txt_data.split("\n")[6:-1]
    for i, line in enumerate(txt_lines):
        print(line)
        items = line.split()
        try:
            energy = float(items[0])
            mu = float(items[1])
            mu_en = float(items[2])
        except:
            try:
                energy = float(items[1]) + 1e-9 # add 1meV to avoid suplicating values
                mu = float(items[2])
                mu_en = float(items[3])
            except:
                print(">>>>>Failed with  Z=%d, line: %s" % (Z, line))
                print(txt_lines)
                raise("Failed with  Z=%d, line: %s" % (Z, line))
                f.close()
        f.write("%g  %g  %g \n" % (energy*1e6, mu, mu_en))

f.close()



    # print(Z, i, energy, mu, mu_en)
