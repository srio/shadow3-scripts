#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
import os
from shutil import copyfile

def getshonecol_from_file(filename,column,nolost=False):
    tmp = Shadow.Beam()
    tmp.load(filename)
    p_a = tmp.getshonecol(column,nolost=nolost)
    return p_a

def read_beam_from_file(filename):
    tmp = Shadow.Beam()
    tmp.load(filename)
    return tmp

def customize_mainframe(oe0,oe1,oe2,oe3,p=30.00,q=10.00,thetaHFM=2.094,thetaVFM=2.094):

    # ;
    # ;  setting angles and distances for a MONTEL system
    # ;
    #
    # ;
    # ; start of customization part ===================================================
    # ;
    #
    # ; define distance p: source to MONTEL (physical), q: MONTEL to detector (physicacal)


    # ; define what are the MONTEL focal distances
    pfoc = p
    qfoc = q
    # ; define the incident angle (grazing, mrad) for MONTEL HFM and VFM, respectively



    # ;
    # ; end of customization part ===================================================
    # ;
    # ; set distances
    oe1.T_SOURCE = p
    oe1.T_IMAGE  = 0.0
    oe2.T_SOURCE = 0.0
    oe2.T_IMAGE  = 0.0
    oe3.T_SOURCE = 0.0
    oe3.T_IMAGE  = q

    # ; set incidence angle
    t1 = 90.0-(thetaVFM*1e-3*180/numpy.pi) # incidence angle in deg for VFM
    t2 = 90.0-(thetaHFM*1e-3*180/numpy.pi) # incidence angle in deg for HFM
    oe1.T_INCIDENCE = t1
    oe2.T_INCIDENCE = t2
    oe3.T_INCIDENCE = 90.0

    oe1.T_REFLECTION = 90.0
    oe2.T_REFLECTION = 90.0
    oe3.T_REFLECTION = 90.0

    oe1.SSOUR = pfoc
    oe2.SSOUR = pfoc
    oe3.SSOUR = pfoc

    oe1.SIMAG = qfoc
    oe2.SIMAG = qfoc
    oe3.SIMAG = qfoc

    oe1.THETA = t1
    oe2.THETA = t2
    oe3.THETA = t1

    return oe0,oe1,oe2,oe3


def get_mainframe():
    #
    # initialize shadow3 source (oe0) and beam
    #

    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 1.0
    oe0.HDIV2 = 1.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.NPOINT = 25000
    oe0.PH1 = 10000.0
    oe0.SIGDIX = 8.84999972e-05
    oe0.SIGDIZ = 7.1999998e-06
    oe0.SIGMAX = 5.7000001e-05
    oe0.SIGMAZ = 1.04e-05
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0

    oe1.DUMMY = 100.0
    oe1.FCYL = 1
    oe1.FHIT_C = 1
    oe1.FMIRR = 2
    oe1.F_DEFAULT = 0
    oe1.RLEN1 = 0.2
    oe1.RLEN2 = 0.2
    oe1.RWIDX1 = 0.02
    oe1.RWIDX2 = 0.02
    oe1.SIMAG = 10.0
    oe1.SSOUR = 30.0
    oe1.THETA = 89.8800226
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 89.8800226
    oe1.T_REFLECTION = 89.8800226
    oe1.T_SOURCE = 30.0

    oe2.ALPHA = 90.0
    oe2.DUMMY = 100.0
    oe2.FCYL = 1
    oe2.FHIT_C = 1
    oe2.FMIRR = 2
    oe2.F_DEFAULT = 0
    oe2.RLEN1 = 0.2
    oe2.RLEN2 = 0.2
    oe2.RWIDX1 = 0.02
    oe2.RWIDX2 = 0.02
    oe2.SIMAG = 10.0
    oe2.SSOUR = 30.0
    oe2.THETA = 89.8800226
    oe2.T_IMAGE = 0.0
    oe2.T_INCIDENCE = 89.8800226
    oe2.T_REFLECTION = 89.8800226
    oe2.T_SOURCE = 0.0

    oe3.ALPHA = 270.0
    oe3.DUMMY = 100.0
    oe3.FCYL = 1
    oe3.FHIT_C = 1
    oe3.FMIRR = 2
    oe3.F_DEFAULT = 0
    oe3.RLEN1 = 0.2
    oe3.RLEN2 = 0.2
    oe3.RWIDX1 = 0.02
    oe3.RWIDX2 = 0.02
    oe3.SIMAG = 10.0
    oe3.SSOUR = 30.0
    oe3.THETA = 89.8800226
    oe3.T_IMAGE = 10.0
    oe3.T_INCIDENCE = 89.8800226
    oe3.T_REFLECTION = 89.8800226
    oe3.T_SOURCE = 0.0


    return oe0,oe1,oe2,oe3

def run_beamline(oe0,oe1,oe2,oe3,iwrite=1):

    beam = Shadow.Beam()
    #Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")


    #
    #run optical element 1
    #
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    beam.traceOE(oe1,1)
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")


    #
    #run optical element 2
    #
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")
    beam.traceOE(oe2,2)
    if iwrite:
        oe2.write("end.02")
        beam.write("star.02")


    #
    #run optical element 3
    #
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")
    beam.traceOE(oe3,3)
    if iwrite:
        oe3.write("end.03")
        beam.write("star.03")


    return beam


if __name__ == "__main__":
    # MACRO 2

    # ; MONTEL run
    # ; implementation of a MONTEL system. It is a KB system on the same piece.
    # ; It is a non-sequential system, because a ray can have its first intercept
    # ; on either the  Vertical focusing mirror (VFM) or the
    # ; Horizontal focusing mirror (HFV) .
    #
    # ;
    # ; To make the simulation, the trick is to define three optical elements VFM,
    # ; HFM, and VFM) on the same position. Shadow runs twice,
    # ; setting one VFM empty: first run uses VFM (oe 1) and HFM (oe2) and
    # ; second run uses HFM (oe 2) and VFM (oe3).
    # ;
    # ; From these two separated runs we select the that first touch the VHF and first touch ; the HFM by looking at their optical path (col 13). The good rays are those
    # ; doubly focuses.
    # ;
    # ; In addition, we compute the rays that suffer a single reflection and no reflections
    # ; at all.
    # ;
    # ; output files:
    # ;     combined.03   contains the rays suffering double reflection (goo ones)
    # ;     combined2.03   contains all rays traced aftyer their physical path, i.e.,
    # ;                    suffering two-bouns, one-bound and zero-bounds. Note that in this
    # ;                    file all rays are "good" (flag>1).
    # ;
    # ; M. Sanchez del Rio, srio@esrf.eu, 2013-07-29 - Python version 2018-03-01
    # ;

    # ; clean output files
    try:
        os.remove("stara.03")
        os.remove("starb.03")
        os.remove("combined.03")
        os.remove("combined2.03")
    except:
        pass


    # ; first run (o.e. 1 & 2)

    oe0, oe1, oe2, oe3 = get_mainframe()
    oe0, oe1, oe2, oe3 = customize_mainframe(oe0, oe1, oe2, oe3)
    oe1.F_REFRAC = 0 #
    oe2.F_REFRAC = 0 #
    oe3.F_REFRAC = 2 # empty

    beam = run_beamline(oe0,oe1,oe2,oe3)

    p_a = getshonecol_from_file("mirr.01",13)
    f_a = getshonecol_from_file('mirr.01',10)

    a = read_beam_from_file('star.03')
    thru_a = numpy.where( ( (f_a < 0) & ( getshonecol_from_file('mirr.02',10) < 0) ) )
    single_a = numpy.where( ( (f_a >= 0) & ( getshonecol_from_file('mirr.02',10) < 0) ) )
    both_a = numpy.where( ( (f_a >= 0) & ( getshonecol_from_file('mirr.02',10) >= 0) ) )

    copyfile("star.03","stara.03")



    # ; second run (o.e. 2 & 3)

    oe0, oe1, oe2, oe3 = get_mainframe()
    oe0, oe1, oe2, oe3 = customize_mainframe(oe0, oe1, oe2, oe3)
    oe1.F_REFRAC = 2 # empty
    oe2.F_REFRAC = 0 #
    oe3.F_REFRAC = 0 #

    beam = run_beamline(oe0,oe1,oe2,oe3)


    p_b = getshonecol_from_file('mirr.02',13)
    f_b = getshonecol_from_file('mirr.02',10)
    b = read_beam_from_file('star.03')
    thru_b   = numpy.where( ( (f_b < 0)  & ( getshonecol_from_file('mirr.03',10) < 0) ) )
    single_b = numpy.where( ( (f_b >= 0) & ( getshonecol_from_file('mirr.03',10) < 0) ) )
    both_b   = numpy.where( ( (f_b >= 0) & ( getshonecol_from_file('mirr.03',10) >= 0) ) )
    copyfile("star.03","starb.03")

    #
    # ; paths and flags
    #
    # ; double-reflection
    c= a.duplicate() # first run to new variable
    # ;copy good rays from second run to new variable
    bound = f_a*0.0 - 1 # define an array to store the bounds: 0=None,1=First,2=Second,3=Both
    for i in range(f_b.size):
        if (f_b[i] > 0):
            bound[i]=2
        if (f_a[i] > 0):
            bound[i]=1
        if ( (f_b[i] > 0) & (p_b[i] < p_a[i])):
            print(i+1,f_a[i],f_b[i],p_a[i],p_b[i],p_b[i]-p_a[i])
            c.rays[i,:] = b.rays[i,:]
            bound[i] = 2
        if (c.rays[i,9] > 0):
            bound[i] = 3
        if ( (f_a[i] < 0) & (f_b[i] < 0)  ):
            bound[i] = 0

    c.write("combined.03")


    #
    # ;
    # ; end of calculations of two-bound rays. Output file: combined.dat
    # ;
    #

    if True:
        # ; if you want to calculate the image given by the rays that suffer only one or
        # ; zero reflections, the following code should be run. It re-runs shadow with
        # ; all oe empty (zero reflections) and setting only one oe. The combined
        # ; results will go to file combined2.03
        #
        #
        #
        #
        # ;recalculates the zero-bound and one-bound rays
        #
        # ; zero bound, all empty
        nbound = numpy.array(numpy.where(bound == 0)).size
        if ( nbound > 1):
            oe0, oe1, oe2, oe3 = get_mainframe()
            oe0, oe1, oe2, oe3 = customize_mainframe(oe0, oe1, oe2, oe3)
            oe1.F_REFRAC = 2
            oe2.F_REFRAC = 2
            oe3.F_REFRAC = 2

            beam = run_beamline(oe0,oe1,oe2,oe3)

            d = read_beam_from_file('star.03')
            c.rays[numpy.where(bound == 0),:] = d.rays[numpy.where(bound == 0),:]


        # ; only intercepts VFM mirror
        nbound = numpy.array(numpy.where(bound == 1)).size
        if ( nbound > 1):
            oe0, oe1, oe2, oe3 = get_mainframe()
            oe0, oe1, oe2, oe3 = customize_mainframe(oe0, oe1, oe2, oe3)
            oe1.F_REFRAC = 0 # empty
            oe2.F_REFRAC = 2 #
            oe3.F_REFRAC = 2 #

            beam = run_beamline(oe0,oe1,oe2,oe3)

            d = read_beam_from_file('star.03')
            c.rays[numpy.where(bound == 1),:] = d.rays[numpy.where(bound == 1),:]



        # ; only intercepts HFM mirror
        nbound = numpy.array(numpy.where(bound == 2)).size
        if ( nbound > 1):
            oe0, oe1, oe2, oe3 = get_mainframe()
            oe0, oe1, oe2, oe3 = customize_mainframe(oe0, oe1, oe2, oe3)
            oe1.F_REFRAC = 2 # empty
            oe2.F_REFRAC = 0 #
            oe3.F_REFRAC = 2 #

            beam = run_beamline(oe0,oe1,oe2,oe3)

            d = read_beam_from_file('star.03')
            c.rays[numpy.where(bound == 2),:] = d.rays[numpy.where(bound == 2),:]

        c.write("combined2.03")


        n0  = numpy.array(numpy.where(bound == 0)).size
        n1v = numpy.array(numpy.where(bound == 1)).size
        n1h = numpy.array(numpy.where(bound == 2)).size
        n2  = numpy.array(numpy.where(bound == 3)).size

        print('rays with no-reflection:           %d '%n0)
        print('rays with single-reflection (VFM): %d '%n1v)
        print('rays with single-reflection (HFM): %d '%n1h)
        print('rays with double-reflection:       %d '%n0)
        print('total:                             %d '%(n0+n1h+n1v+n2))

        #
        # display results
        #
        x,z = read_beam_from_file("combined2.03").getshcol([1,3])
        from srxraylib.plot.gol import plot_scatter
        plot_scatter(x,z,title='tracing all rays',show=False)


    x,z = read_beam_from_file("combined.03").getshcol([1,3])
    from srxraylib.plot.gol import plot_scatter
    plot_scatter(x,z,title='rays with two reflections')



