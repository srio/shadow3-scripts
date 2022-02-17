import numpy
import scipy.constants as codata



def check_ID19_W150b(show=1):

    # 19,W150b,Middle,150.0,1030.00, 26.500, 2.6479, 0.0000, 0.8699, 1.2400, 1.0000, 1.1300, ,
    A = [2.6479, 0.0000, 0.8699, 1.2400, 1.0000, 1.1300]
    id_period_mm = 150.0
    min_gap = 26.500


    gap_mm = numpy.linspace(min_gap, 2*min_gap, 100)


    Bmax = numpy.zeros_like(gap_mm)

    Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * 1 * gap_mm / id_period_mm)
    Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * 2 * gap_mm / id_period_mm)
    Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * 3 * gap_mm / id_period_mm)



    K = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)
    K2 = (93.4 * 0.01 * (id_period_mm/10) * (
                2.3333 * numpy.exp(-0.02473 * gap_mm) + 1.189 * numpy.exp(-0.059691 * gap_mm)))

    print(Bmax, K)

    from srxraylib.plot.gol import plot
    plot(gap_mm, K, gap_mm, K2, xtitle="gap [mm]", ytitle="K", legend=['now','gwen'],title="ID19 W150b", show=show)

def check_ID06_CPMU(show=1):
    # 6,CPMU18a,UpStream,18.0,2000.00,  6.000, 2.6600, 0.0000, 0.8720, 1.0826, 1.0000, 1.1456, ,
    A = [2.6600, 0.0, 0.8720, 1.0826, 1.0000, 1.1456]
    id_period_mm = 18.0
    min_gap = 5.0


    gap_mm = numpy.linspace(min_gap, 2*min_gap, 100)


    Bmax = numpy.zeros_like(gap_mm)
    Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * gap_mm / id_period_mm)
    Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * gap_mm / id_period_mm)
    Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * gap_mm / id_period_mm)
    K = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)

    Bmax = numpy.zeros_like(gap_mm)
    Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * 1 * gap_mm / id_period_mm)
    Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * 2 * gap_mm / id_period_mm)
    Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * 3 * gap_mm / id_period_mm)
    K2 = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)


    Bmax3 = 3.29 * numpy.exp(-numpy.pi * 1.095 * gap_mm / id_period_mm)
    K3 = Bmax3 * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)


    from srxraylib.plot.gol import plot
    plot(gap_mm, K, gap_mm, K2, gap_mm, K3, xtitle="gap [mm]", ytitle="K", legend=['old', 'new','chavanne generic'],
         title="ID06 CPMU", show=show)


def check_ID11_CMPU(show=1):

    A = [2.978, 0.0, 0.988, 1.0826, 1.0000, 1.1456]
    id_period_mm = 18.0
    min_gap = 4.0


    gap_mm = numpy.linspace(min_gap, 2*min_gap, 100)


    Bmax = numpy.zeros_like(gap_mm)
    Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * gap_mm / id_period_mm)
    Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * gap_mm / id_period_mm)
    Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * gap_mm / id_period_mm)
    K = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)


    Bmax = numpy.zeros_like(gap_mm)
    Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * 1 * gap_mm / id_period_mm)
    Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * 2 * gap_mm / id_period_mm)
    Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * 3 * gap_mm / id_period_mm)
    K2 = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)


    Bmax3 = 3.29 * numpy.exp(-numpy.pi * 1.095 * gap_mm / id_period_mm)
    K3 = Bmax3 * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)

    from srxraylib.plot.gol import plot
    plot(gap_mm, K, gap_mm, K2, gap_mm, K3, xtitle="gap [mm]", ytitle="K", legend=['old', 'new','chavanne generic'],
         title="ID11 CPMU", show=show)

def check_ID15_CMPU(show=1):

    A = [3.2077, 0.0, 0.9944, 1.0771, 1.0000, 1.1184]
    id_period_mm = 18.0
    min_gap = 4.0


    gap_mm = numpy.linspace(min_gap, 2*min_gap, 100)


    Bmax = numpy.zeros_like(gap_mm)
    Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * gap_mm / id_period_mm)
    Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * gap_mm / id_period_mm)
    Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * gap_mm / id_period_mm)
    K = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)


    Bmax = numpy.zeros_like(gap_mm)
    Bmax += A[0] * numpy.exp(-numpy.pi * A[3] * 1 * gap_mm / id_period_mm)
    Bmax += A[1] * numpy.exp(-numpy.pi * A[4] * 2 * gap_mm / id_period_mm)
    Bmax += A[2] * numpy.exp(-numpy.pi * A[5] * 3 * gap_mm / id_period_mm)
    K2 = Bmax * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)


    Bmax3 = 3.29 * numpy.exp(-numpy.pi * 1.095 * gap_mm / id_period_mm)
    K3 = Bmax3 * (id_period_mm * 1e-3) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)

    from srxraylib.plot.gol import plot
    plot(gap_mm, K, gap_mm, K2, gap_mm, K3, xtitle="gap [mm]", ytitle="K", legend=['old', 'new','chavanne generic'],
         title="ID15 CPMU", show=show)



if __name__ == "__main__":
    check_ID19_W150b(show=0)
    check_ID06_CPMU( show=0)
    check_ID11_CMPU( show=0)
    check_ID15_CMPU( show=1)