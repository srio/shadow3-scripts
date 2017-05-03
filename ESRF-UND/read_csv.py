
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

        if year == "2017":
            src1._energy_in_GeV = 6.037

            if numpy.mod(int(Section),2):
                source = "LowBeta"
                src1.set_sigmas_all(sigma_x=0.0,sigma_xp=0.0,sigma_y=0.0,sigma_yp=0.0)
                print("Section: %s Lb"%Section)
            else:
                source = "HighBeta"
                src1.set_sigmas_all(sigma_x=0.0,sigma_xp=0.0,sigma_y=0.0,sigma_yp=0.0)
                print("Section: %s Hb"%Section)
        else:
            source = "EBS"
            src1._energy_in_GeV = 6.0
            src1.set_sigmas_all(sigma_x=0.0,sigma_xp=0.0,sigma_y=0.0,sigma_yp=0.0)
            print("Section: %s EBS"%Section)



        UlistS = Ulist.split("/")
        print(">>>",(UlistS[2]))
        #
        src = LightSource("ESRF_%s_ID%s_%s_%s_%s"%(source,Section,UlistS[0],UlistS[1],UlistS[2]),src1,src2)

        print(src.info())


