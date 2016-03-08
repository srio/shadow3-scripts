#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
from srxraylib.sources import srfunc
import matplotlib
from matplotlib import pylab as plt
import os

matplotlib.rcParams.update({'font.size': 8})


def wiggler_preprocessor(ener_gev=6.0,e_min=5000.0,e_max=5005.0,file_field="",plot_trajectories=0,
                        shift_x_flag=0,shift_x_value=0.0,shift_betax_flag=0,shift_betax_value=0.0):

    (traj, pars) = srfunc.wiggler_trajectory(b_from=1,
                                             inData=file_field,
                                             nPer=1,
                                             nTrajPoints=501,
                                             ener_gev=ener_gev,
                                             per=None,
                                             kValue=None,
                                             trajFile="tmp.traj",
                                             shift_x_flag=shift_x_flag,
                                             shift_x_value=shift_x_value,
                                             shift_betax_flag=shift_betax_flag,
                                             shift_betax_value=shift_betax_value)


    data = numpy.loadtxt("tmp.traj",skiprows=15)

    fig = plt.figure(1)

    fig.add_subplot(221)
    plt.plot(data[:,1],data[:,7])
    plt.title("Magnetic Field "+file_field)
    plt.xlabel("Y [m]")
    plt.ylabel("B [T]")

    fig.add_subplot(222)
    plt.plot(data[:,1],data[:,0])
    plt.title("Electron trajectory")
    plt.xlabel("Y [m]")
    plt.ylabel("X [m]")

    fig.add_subplot(223)
    plt.plot(data[:,1],data[:,3])
    plt.title("Electron velocity")
    plt.xlabel("Y [m]")
    plt.ylabel("betaX")

    fig.add_subplot(224)
    plt.plot(data[:,1],data[:,6])
    plt.title("Electron curvature")
    plt.xlabel("Y [m]")
    plt.ylabel("curvature [m^-1]")



    if plot_trajectories:
        plt.show()

    fig.savefig('sw_'+file_field+'.png')
    plt.close(fig)

    srfunc.wiggler_cdf(traj,
                       enerMin=e_min,
                       enerMax=e_max,
                       enerPoints=1001,
                       outFile="xshwig.sha",
                       elliptical=False)

def wiggler_source(ener_gev=6.0,e_min=5000.0,e_max=5005.0,iwrite=1,emittance=0):

    beam = Shadow.Beam()
    oe0 = Shadow.Source()

    oe0.BENER = 6.0
    oe0.CONV_FACT = 100.0
    oe0.FDISTR = 0
    oe0.FILE_TRAJ = b'xshwig.sha'
    oe0.FSOUR = 0
    oe0.FSOURCE_DEPTH = 0
    oe0.F_COLOR = 0
    oe0.F_PHOT = 0
    oe0.F_WIGGLER = 1
    oe0.HDIV1 = 1.0
    oe0.HDIV2 = 1.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.NPOINT = 20000
    oe0.NTOTALPOINT = 0
    oe0.N_COLOR = 0
    oe0.PH1 = e_min
    oe0.PH2 = e_max
    oe0.POL_DEG = 0.0
    oe0.SIGMAY = 0.0
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0



    if emittance == 0:
        oe0.SIGMAX = 0.0
        oe0.SIGMAZ = 0.0
        oe0.EPSI_X = 0.0
        oe0.EPSI_Z = 0.0
        oe0.EPSI_DX = 0.0
        oe0.EPSI_DZ = 0.0
    else:
        # oe0.SIGMAX = 0.001373
        # oe0.SIGMAZ = 3.64e-04
        # oe0.EPSI_X =  oe0.SIGMAX * 16.524e-6 # 2.28e-08
        # oe0.EPSI_Z =  oe0.SIGMAZ * 1.623e-6  # 5e-10
        # oe0.EPSI_DX = 65.1
        # oe0.EPSI_DZ = -28.55
        oe0.SIGMAX = 0.0008757
        oe0.SIGMAZ = 0.0001647
        oe0.EPSI_X =  oe0.SIGMAX * 24.72e-6 # 2.16e-08
        oe0.EPSI_Z =  oe0.SIGMAZ * 3.036e-6  # 5e-10
        oe0.EPSI_DX = 89.4
        oe0.EPSI_DZ = -104.8



    #Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    return beam,oe0

def bm_source(ener_gev=6.0,e_min=5000.0,e_max=5005.0,iwrite=0,emittance=0):

    beam = Shadow.Beam()
    oe0 = Shadow.Source()



    oe0.BENER = 6.04
    if emittance == 0:
        oe0.SIGMAX = 0.0
        oe0.SIGMAZ = 0.0
        oe0.EPSI_X = 0.0
        oe0.EPSI_Z = 0.0
        oe0.EPSI_DX = 0.0
        oe0.EPSI_DZ = 0.0

    else:
        oe0.SIGMAX = 77.9e-4
        oe0.SIGMAZ = 12.9e-4
        oe0.EPSI_X = oe0.SIGMAX * 110.9e-6  # 8.6e-07
        oe0.EPSI_Z = oe0.SIGMAZ * 0.5e-6    # 6.45e-10
        oe0.EPSI_DX = 0.0
        oe0.EPSI_DZ = 0.0

    oe0.FDISTR = 4
    oe0.FSOURCE_DEPTH = 4
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.001
    oe0.HDIV2 = 0.001
    oe0.NCOL = 0
    oe0.N_COLOR = 0
    oe0.PH1 = e_min
    oe0.PH2 = e_max
    oe0.POL_DEG = 0.0
    oe0.R_ALADDIN = 2353
    oe0.R_MAGNET = 23.53
    oe0.SIGDIX = 0.0
    oe0.SIGDIZ = 0.0

    oe0.NPOINT = 20000
    oe0.SIGMAY = 0.0
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0


    #Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    return beam,oe0

def focusing_mirror(beam,foc="1:1",grazing_theta_mrad=3.0,iwrite=0):


    if foc == "1:1":
        dist_p = 3000.0
        dist_q = 3000.0
        fmirr = 3  # toroid
    elif foc == "3:1":
        dist_p = 4500.0
        dist_q = 1500.0
        fmirr = 2 # ellipsoid
    else:
        raise "Problems..."

    oe1 = Shadow.OE()
    oe1.DUMMY = 1.0
    oe1.FHIT_C = 1
    oe1.FMIRR = fmirr
    oe1.FWRITE = 3
    oe1.RLEN1 = 50.0
    oe1.RLEN2 = 50.0
    oe1.RWIDX1 = 10.0
    oe1.RWIDX2 = 10.0
    oe1.T_IMAGE = dist_q
    oe1.T_INCIDENCE = 90.0 - grazing_theta_mrad * 1e-3 * 180 / numpy.pi
    oe1.T_REFLECTION = oe1.T_INCIDENCE
    oe1.T_SOURCE = dist_p


    # run optical element 1

    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    beam.traceOE(oe1,1)
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")
    return beam,oe1


def focusing_ideal(beam,iwrite=0):
    # implements an infinite elliptical mirror at 45 deg.
    oe1 = Shadow.OE()
    oe1.DUMMY = 1.0
    oe1.FMIRR = 2
    oe1.FWRITE = 3
    oe1.T_IMAGE = 3000.0
    oe1.T_INCIDENCE = 45.0
    oe1.T_REFLECTION = 45.0
    oe1.T_SOURCE = 3000.0

    # run optical element 1

    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    beam.traceOE(oe1,1)
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    return beam,oe1

if __name__ == "__main__":

    wigglers =  ["2PAcut"]       # ["BM","3Pcut","2PAcut","2PBcut","1Pcut","3P","2PA","2PB","1P","0P"] #
    energies =  [5]              #  [5,10,20,40,80]                                                    #
    emittance = 1

    # initialize
    plot_trajectories = 0
    txt_all = ""
    txt_latex = ""
    plot_index = -1

    xrange_1 = [-0.015,0.015]
    yrange_1 = [-0.002,0.002]
    xrange_3 = xrange_1 # [-0.015,0.015]
    yrange_3 = yrange_1 # [-0.01,0.01]

    for root_file in wigglers:
        if emittance:
            title0 = "Finite Emittance; magnetic field: %s"%(root_file)
        else:
            title0 = "Zero Emittance; magnetic field: %s"%(root_file)
        txt_all += "\n == %s ==\n\n"%(title0)

        file_field = "SW_"+root_file+".txt"
        #plot trajectories
        txt_all += "\n\n [%s] \n\n"%('sw_'+file_field+'.png')


        for energy_kev in energies:
            # start ENERGY loop HERE #######################


            e_min = energy_kev * 1000.0
            e_max = e_min + e_min/1000.0 # 0.1%bw

            if emittance:
                title = "Finite Emittance; magnetic field: %s; E=%d keV"%(root_file,energy_kev)
                png_label="_withemitt_"
            else:
                title = "Zero Emittance; magnetic field: %s; E=%d keV"%(root_file,energy_kev)
                png_label="_zeroemitt_"


            txt_all += "\n === %s ===\n\n"%(title)


            #
            # source
            #

            # preprocessor ["BM","3Pcut","2PAcut","2PBcut","1Pcut"]
            if root_file == "BM":
                pass
            #
            # first, only the short IDs
            #
            elif root_file == "1Pcut":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=4,shift_x_flag=4)
            elif root_file == "2PAcut":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=1,shift_x_flag=1)
            elif root_file == "2PBcut":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=1,shift_x_flag=1)
            elif root_file == "3Pcut":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=0,shift_x_flag=1)
            #
            # now, also with side BMs
            #
            elif root_file == "0P":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=4,shift_x_flag=4)
            elif root_file == "1P":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=4,shift_x_flag=4)
            elif root_file == "2PA":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=4,shift_x_flag=4)

            elif root_file == "2PB":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=4,shift_x_flag=4)
            elif root_file == "3P":
                wiggler_preprocessor(e_min=e_min,e_max=e_max,file_field=file_field,
                                     plot_trajectories=plot_trajectories,shift_betax_flag=4,shift_x_flag=4)
            else:
                raise "Error: root_file not valid"

            #shadow source
            if root_file == "BM":
                beam, oe0 = bm_source(e_min=e_min,e_max=e_max,emittance=emittance)
            else:
                beam, oe0 = wiggler_source(e_min=e_min,e_max=e_max,emittance=emittance)

            tmp0 = oe0.sourcinfo()
            # plot divergences
            subtitle = "Divergence space"
            tkt0 = Shadow.ShadowTools.plotxy_gnuplot(beam,4,6,xrange=[-0.003,0.003],yrange=[-250e-6,250e-6],nolost=1,
                    ref="Yes",nbins=201,ps=1,viewer="ls ",title="Div space; Units: rad; "+subtitle,
                    calculate_widths=2)
            if tkt0["fw25%m_h"] == None: tkt0["fw25%m_h"] = -1e-4
            if tkt0["fw25%m_v"] == None: tkt0["fw25%m_v"] = -1e-4
            if tkt0["fwhm_h"] == None:   tkt0["fwhm_h"] = -1e-4
            if tkt0["fwhm_v"] == None:   tkt0["fwhm_v"] = -1e-4
            if tkt0["fw75%m_h"] == None: tkt0["fw75%m_h"] = -1e-4
            if tkt0["fw75%m_v"] == None: tkt0["fw75%m_v"] = -1e-4

            txt_all += "\n ==== %s ====\n\n"%(subtitle)
            txt_all += "\n\n"
            txt_all += " | FW25M  | H: %6.3f urad | V: %6.3f urad \n"%(1e6*tkt0["fw25%m_h"],1e6*tkt0["fw25%m_v"])
            txt_all += " | FWHM   | H: %6.3f urad | V: %6.3f urad \n"%(1e6*tkt0["fwhm_h"],  1e6*tkt0["fwhm_v"])
            txt_all += " | FW75M  | H: %6.3f urad | V: %6.3f urad \n"%(1e6*tkt0["fw75%m_h"],1e6*tkt0["fw75%m_v"])

            plot_index += 1
            file_png = "sw%s%000d.png"%(png_label,plot_index)
            os.system("convert -rotate 90 plotxy.ps %s"%(file_png))
            txt_all += "\n\n [%s] \n\n"%(file_png)

            # plot top view
            subtitle = "Emission intensity - electron trajectory (top view)"
            txt_all += "\n ==== %s ====\n\n"%(subtitle)
            txt_all += "\n\n"
            Shadow.ShadowTools.plotxy_gnuplot(beam,2,1,nbins=201,ps=1,viewer="ls ",title="Units=cm; %s"%subtitle)


            plot_index += 1
            file_png = "sw%s%000d.png"%(png_label,plot_index)
            os.system("convert -rotate 90 plotxy.ps %s"%(file_png))
            txt_all += "\n\n [%s] \n\n"%(file_png)


            #
            # ideal focusing
            #

            beam,oe1 = focusing_ideal(beam)

            subtitle = "Ideal 1:1 focusing"

            # plot divergences
            # xrange=[-500e-4/4,500e-4/4],yrange=[-500e-4/16,500e-4/16]
            tkt1 = Shadow.ShadowTools.plotxy_gnuplot(beam,1,3,xrange=xrange_1, yrange=yrange_1,
                                                     nolost=1,ref="Yes",nbins=201,ps=1,viewer="ls ",
                                                     title="Real space; Units=cm; %s"%subtitle,
                                                     calculate_widths=2)
            if tkt1["fw25%m_h"] == None: tkt1["fw25%m_h"] = -1e-4
            if tkt1["fw25%m_v"] == None: tkt1["fw25%m_v"] = -1e-4
            if tkt1["fwhm_h"] == None:   tkt1["fwhm_h"] = -1e-4
            if tkt1["fwhm_v"] == None:   tkt1["fwhm_v"] = -1e-4
            if tkt1["fw75%m_h"] == None: tkt1["fw75%m_h"] = -1e-4
            if tkt1["fw75%m_v"] == None: tkt1["fw75%m_v"] = -1e-4


            txt_all += "\n ==== %s ====\n\n"%(subtitle)
            txt_all += "\n\n"
            txt_all += " | FW25M  | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt1["fw25%m_h"],1e4*tkt1["fw25%m_v"])
            txt_all += " | FWHM   | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt1["fwhm_h"],  1e4*tkt1["fwhm_v"])
            txt_all += " | FW75M  | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt1["fw75%m_h"],1e4*tkt1["fw75%m_v"])

            plot_index += 1
            file_png = "sw%s%000d.png"%(png_label,plot_index)
            os.system("convert -rotate 90 plotxy.ps %s"%(file_png))
            txt_all += "\n\n [%s] \n\n"%(file_png)


            #
            # 1:1 focusing
            #
            beam = None
            if root_file == "BM":
                beam, oe0 = bm_source(e_min=e_min,e_max=e_max,emittance=emittance)
            else:
                beam, oe0 = wiggler_source(e_min=e_min,e_max=e_max,emittance=emittance)


            beam,oe1 = focusing_mirror(beam,foc="1:1",grazing_theta_mrad=3.0)

            subtitle = "Toroid 1:1 focusing"

            tkt2 = Shadow.ShadowTools.plotxy_gnuplot(beam,1,3,xrange=xrange_1, yrange=yrange_1,
                                                     nolost=1,ref="Yes",nbins=301,ps=1,viewer="ls ",
                                                     title="Real space; Units=cm; %s"%subtitle,
                                                     calculate_widths=2)
            if tkt2["fw25%m_h"] == None: tkt2["fw25%m_h"] = -1e-4
            if tkt2["fw25%m_v"] == None: tkt2["fw25%m_v"] = -1e-4
            if tkt2["fwhm_h"]   == None: tkt2["fwhm_h"] = -1e-4
            if tkt2["fwhm_v"]   == None: tkt2["fwhm_v"] = -1e-4
            if tkt2["fw75%m_h"] == None: tkt2["fw75%m_h"] = -1e-4
            if tkt2["fw75%m_v"] == None: tkt2["fw75%m_v"] = -1e-4

            txt_all += "\n ==== %s ====\n\n"%(subtitle)
            txt_all += "\n\n"
            txt_all += " | FW25M  | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt2["fw25%m_h"],1e4*tkt2["fw25%m_v"])
            txt_all += " | FWHM   | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt2["fwhm_h"],  1e4*tkt2["fwhm_v"])
            txt_all += " | FW75M  | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt2["fw75%m_h"],1e4*tkt2["fw75%m_v"])

            plot_index += 1
            file_png = "sw%s%000d.png"%(png_label,plot_index)
            os.system("convert -rotate 90 plotxy.ps %s"%(file_png))
            txt_all += "\n\n [%s] \n\n"%(file_png)
            #

            #
            # 3:1 focusing
            #
            beam = None
            if root_file == "BM":
                beam, oe0 = bm_source(e_min=e_min,e_max=e_max,emittance=emittance)
            else:
                beam, oe0 = wiggler_source(e_min=e_min,e_max=e_max,emittance=emittance)
            beam,oe1 = focusing_mirror(beam,foc="3:1",grazing_theta_mrad=3.0)

            subtitle = "Ellipsoid 3:1 focusing"
            # xrange=[-500e-4/2,500e-4/2],yrange=[-500e-4/4,500e-4/4]
            tkt3 = Shadow.ShadowTools.plotxy_gnuplot(beam,1,3, xrange=xrange_3, yrange=yrange_3,
                                                     nolost=1,ref="Yes",nbins=151,ps=1,viewer="ls ",
                                                     title="Real space; Units=cm; %s"%subtitle,
                                                     calculate_widths=2)

            if tkt3["fw25%m_h"] == None: tkt3["fw25%m_h"] = -1e-4
            if tkt3["fw25%m_v"] == None: tkt3["fw25%m_v"] = -1e-4
            if tkt3["fwhm_h"] == None:   tkt3["fwhm_h"] = -1e-4
            if tkt3["fwhm_v"] == None:   tkt3["fwhm_v"] = -1e-4
            if tkt3["fw75%m_h"] == None: tkt3["fw75%m_h"] = -1e-4
            if tkt3["fw75%m_v"] == None: tkt3["fw75%m_v"] = -1e-4

            txt_all += "\n ==== %s ====\n\n"%(subtitle)
            txt_all += "\n\n"
            txt_all += " | FW25M  | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt3["fw25%m_h"],1e4*tkt3["fw25%m_v"])
            txt_all += " | FWHM   | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt3["fwhm_h"],  1e4*tkt3["fwhm_v"])
            txt_all += " | FW75M  | H: %6.3f um | V: %6.3f um \n"%(1e4*tkt3["fw75%m_h"],1e4*tkt3["fw75%m_v"])

            plot_index += 1
            file_png = "sw%s%000d.png"%(png_label,plot_index)
            os.system("convert -rotate 90 plotxy.ps %s"%(file_png))
            txt_all += "\n\n [%s] \n\n"%(file_png)

            if root_file != "0P":
                txt_latex += "%d  &  %s  & %d  & %d & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f %s%s \n"%(energy_kev,root_file,\
                                    1e6*tkt0["fwhm_h"],1e6*tkt0["fwhm_v"],\
                                    1e4*tkt1["fwhm_h"],1e4*tkt1["fwhm_v"],\
                                    1e4*tkt2["fwhm_h"],1e4*tkt2["fwhm_v"],\
                                    1e4*tkt3["fwhm_h"],1e4*tkt3["fwhm_v"],\
                                    "\\","\\")
        txt_latex += "\hline \n"

    #
    # dump results
    #
    if emittance:
        file_latex = "table2.txt"
    else:
        file_latex = "table1.txt"

    f = open(file_latex,"w")
    f.write(txt_latex)
    f.close()
    print("File written to disk: %s"%(file_latex))


    if emittance:
        root_file_out = 'short_wigglers_with_emittance'
    else:
        root_file_out = 'short_wigglers_no_emittance'

    file_out = root_file_out+'.t2t'
    f = open(file_out,"w")
    f.write(txt_all)
    f.close()
    print("File written to disk: ",file_out)

    os.system("./txt2tags -t html --toc --enum-title %s"%file_out)
    print("All done.")
    print(tmp0)

    # if plot_trajectories: plt.show()
    #os.system("open %s.html"%root_file_out)



