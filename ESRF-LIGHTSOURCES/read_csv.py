


from syned.storage_ring.electron_beam import ElectronBeam
from syned.storage_ring.magnetic_structures.undulator import Undulator
from syned.storage_ring.light_source import LightSource
import numpy




f = open('ID_data_jan_2017_srio.csv','r')
txt = f.read().split('\n')


#
for i,line in enumerate(txt):
    if i != 0:

        Ulist,Carriage,Section,Position,Period,Length,Kzmax,Kxmax,year = line.split(",")


        src1 = ElectronBeam()
        src2 = Undulator()


        src1._current = 0.2
        src1._number_of_bunches = 1

        if year == "2017":
            src1._energy_spread = 0.0011
            src1._energy_in_GeV = 6.037

            if numpy.mod(int(Section),2):
                source = "LowBeta"
                print("Section: %s Lb"%Section)
                # src1.set_sigmas_all(37.4e-6,106.9e-6,3.5e-6,1.2e-6)
                src1.set_sigmas_all(5.03973e-05,0.000107192,3.44459e-06,1.16124e-06)
            else:
                source = "HighBeta"
                print("Section: %s Hb"%Section)
                # src1.set_sigmas_all(387.8e-6,10.3e-6,3.5e-6,1.2e-6)
                src1.set_sigmas_all(0.000414971,1.03149e-05,3.43353e-06,1.16498e-06)

        else:
            src1._energy_spread = 9.3339e-04
            # src1._energy_spread = 9e-04
            src1._energy_in_GeV = 6.0
            source = "EBS"
            print("Section: %s EBS"%Section)
            # src1.set_sigmas_all(27.2e-6,5.2e-6,3.4e-6,1.4e-6)
            src1.set_sigmas_all(3.01836e-05,4.36821e-06,3.63641e-06,1.37498e-06)

        src2._K_horizontal = float(Kxmax)
        src2._K_vertical = float(Kzmax)
        src2._period_length = float(Period)*1e-3
        src2._number_of_periods = int( float(Length) / src2._period_length)

        UlistS = Ulist.split("/")
        print(">>>",(UlistS[2]))
        #
        name = "ESRF_ID%s_%s_%s_%s"%(Section,source,UlistS[1],UlistS[2])
        src = LightSource(name,src1,src2)

        print(src.info())
        src.to_json(name+'.json')


