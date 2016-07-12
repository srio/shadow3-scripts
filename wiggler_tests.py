import numpy
import matplotlib.pylab as plt
from scipy import interpolate
import Shadow

def read_file_cdf(file_cdf="xshwig.sha"):

    f = open(file_cdf,"r")
    firstline = f.readline()
    f.close()

    fl = firstline.split(" ")[0:-1]

    print(len(firstline),fl)
    np           = int(fl[0])
    step         = float(fl[1])
    bener        = float(fl[2])
    r_min        = float(fl[3])
    r_max        = float(fl[4])
    energy_start = float(fl[5])
    energy_end   = float(fl[6])

    # print("np           = %d"%( np ))
    # print("step         = %f"%( step ))
    # print("bener        = %f"%( bener ))
    # print("rad_min      = %f"%( r_max ))
    # print("rad_min      = %f"%( r_min ))
    # print("energy_start = %f"%( energy_start ))
    # print("energy_end   = %f"%( energy_end ))

    header = {}
    header["NP_TRAJ"] = np
    header["PATH_STEP"] = step
    header["BENER"] = bener
    header["RAD_MIN"] = r_max
    header["RAD_MAX"] = r_min
    header["PH1"] = energy_start
    header["PH2"] = energy_end


    alldata = numpy.loadtxt(file_cdf,skiprows=1)

    return header,alldata


def plot_cdf(alldata,label="cdf"):

    x     = alldata[:,0].reshape(-1)
    y     = alldata[:,1].reshape(-1)
    cdf   = alldata[:,2].reshape(-1)
    angle = alldata[:,3].reshape(-1)
    curv  = alldata[:,4].reshape(-1)

    if label == "x":
        ordinates_to_plot = x
    elif label == "y":
        ordinates_to_plot = y
    elif label == "cdf":
        ordinates_to_plot = cdf
    elif label == "angle":
        ordinates_to_plot = angle
    elif label == "curv":
        ordinates_to_plot = curv

    plt.figure(1)

    plt.subplot(221)
    plt.plot(y,1e6*x,'b',label=label)
    plt.title("Trajectory")
    plt.xlabel("y [m]")
    plt.ylabel("x [$\mu m$]")


    plt.subplot(222)
    plt.plot(y,cdf,'b',label=label)
    plt.title("CDF")
    plt.xlabel("y [m]")
    plt.ylabel("cdf")

    plt.subplot(223)
    plt.plot(y,angle,'b',label=label)
    plt.title("Angle")
    plt.xlabel("y [m]")
    plt.ylabel("angle [rad]")

    plt.subplot(224)
    plt.plot(y,curv,'b',label=label)
    plt.title("Curvature")
    plt.xlabel("y [m]")
    plt.ylabel("curvature [m^-1]")

    plt.show()

def plot_xy(x,y,label="y(x)",symbol='b'):


    plt.figure(1)

    plt.subplot(111)
    plt.plot(x,y,symbol,label=label)
    plt.title(label)
    plt.xlabel("")
    plt.ylabel("")

    plt.show()

if __name__ == "__main__":

    (cdf_header,cdf_data) = read_file_cdf(file_cdf="xshwig.sha")

    print(cdf_data.shape)

    # plot_cdf(cdf_data ,label="angle")

    NP_TRAJ    = cdf_header["NP_TRAJ"]
    PATH_STEP  = cdf_header["PATH_STEP"]
    BENER      = cdf_header["BENER"]
    RAD_MIN    = cdf_header["RAD_MIN"]
    RAD_MAX    = cdf_header["RAD_MAX"]
    PH1        = cdf_header["PH1"]
    PH2        = cdf_header["PH2"]

    print(NP_TRAJ,PATH_STEP,BENER,RAD_MIN,RAD_MAX,PH1,PH2)

    CONV_FACT = 100.0

    cdf_path = numpy.linspace(start=0.0,stop=NP_TRAJ,num=NP_TRAJ) * PATH_STEP



    print(cdf_path.shape)


    # seed_y = numpy.zeros( (5,NP_TRAJ) ) # cdf
    # y_x    = numpy.zeros( (5,NP_TRAJ) ) #
    # y_xpri = numpy.zeros( (5,NP_TRAJ) ) #
    # y_curv = numpy.zeros( (5,NP_TRAJ) ) #
    # y_path = numpy.zeros( (5,NP_TRAJ) ) # trajectory
    # # y_z    = numpy.zeros( (5,NP_TRAJ) )
    # # y_zpri = numpy.zeros( (5,NP_TRAJ) )




    y_x    = interpolate.interp1d(cdf_data[:,1]*CONV_FACT,cdf_data[:,0]*CONV_FACT)
    y_xpri = interpolate.interp1d(cdf_data[:,1]*CONV_FACT,cdf_data[:,3])
    y_curv = interpolate.interp1d(cdf_data[:,1]*CONV_FACT,cdf_data[:,4])
    y_path = interpolate.interp1d(cdf_data[:,1]*CONV_FACT,cdf_path*CONV_FACT)
    # seed_y = interpolate.InterpolatedUnivariateSpline(cdf_data[:,2],cdf_data[:,1]*CONV_FACT)
    seed_y = interpolate.interp1d(cdf_data[:,2],cdf_data[:,1]*CONV_FACT)

    # seed_y = interpolate.InterpolatedUnivariateSpline(cdf_data[:,1]*CONV_FACT,cdf_data[:,2])
    # print(cdf_data[:,1],cdf_path)



    src = Shadow.Source()
    src.load("start.00")

    TORAD = numpy.pi / 180.0

    src.F_WIGGLER = 1
    src.F_PHOT  = 0
    src.F_COLOR  = 3
    src.FSOUR  = 3
    src.FDISTR  = 4

    if src.FDISTR == 4:
        src.F_COHER = 0
        if src.R_ALADDIN <= 0.0:
            src.POL_ANGLE = -90.0
        else:
            src.POL_ANGLE = 90.0

    src.POL_ANGLE *= TORAD

    if src.FSOUR == 3:
        EPSI_XOLD = src.EPSI_X
        EPSI_ZOLD = src.EPSI_Z
        if src.SIGMAX != 0.0:
            EPSI_X =   src.EPSI_X / src.SIGMAX
        else:
            EPSI_X = 0.0

        if src.SIGMAZ != 0.0:
            EPSI_Z = src.EPSI_Z / src.SIGMAZ
        else:
            EPSI_Z = 0.0

    GRID = numpy.random.rand(6,src.NPOINT)
    GRID[4,:] = 0.0

    print("GRID shape: ",GRID.shape)


    ARG_Y = GRID[2,:] # numpy.linspace(0.0,1,1000) # GRID[2,:]

    YTRAJ = seed_y(ARG_Y)


    print(ARG_Y,YTRAJ)

    plot_xy(ARG_Y,YTRAJ,symbol="o")



