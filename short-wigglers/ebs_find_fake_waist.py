"""
This script is used to find the "fake" waists for the ESRF-EBS lattice around the short wigglers positions

Input: file  /scisoft/users/srio/Working/ATCOLLAB/esrf2_wiggler.spec  created using AT


"""

from PyMca5.PyMcaIO import specfile
import matplotlib.pylab as plt
import numpy

def sign_change(a):
    # http://stackoverflow.com/questions/2652368/how-to-detect-a-sign-change-for-elements-in-a-numpy-array
    # a = array([1,1,-1,-2,-3,4,5])
    asign = numpy.sign(a)
    sz = asign == 0
    while sz.any():
        asign[sz] = numpy.roll(asign, 1)[sz]
        sz = asign == 0
    signchange = ((numpy.roll(asign, 1) - asign) != 0).astype(int)
    # print signchange
    # array([0, 0, 1, 0, 0, 1, 0])
    signchange[0] = 0
    return signchange

def straight_line_coefficients(x0,y0,x1,y1):
    # returns (m,b) of y=m x + b ; works for arrays!
    m = (y1 - y0) / (x1 - x0)
    b = y1 - m * x1
    return (m,b)

def get_roots(a):
    #returns an array with the indices (float) where the zeroes are found (linear interpolation used)
    iroots = numpy.array(numpy.where(sign_change(a) == 1)[0])
    m, b = straight_line_coefficients(iroots-1,a[iroots-1],iroots,a[iroots])
    return -b/m

def interpolate_at_indices(x,indices=[1.1,2.2]):
    i1 = (numpy.floor(indices)).astype(int)
    i2 =  (numpy.ceil(indices)).astype(int)
    # print(">><<",index,i1,i2)
    val = x[i1] + (indices - i1) * (x[i2] - x[i1])
    return val



class at_data(object):

    def __init__(self): #,epsilon_x=134e-12,epsilon_y=5e-12,delta=0.001,at_file=""):
        self.at_file = ""
        self.epsilon_x = 0.0
        self.epsilon_y = 0.0
        self.delta = 0.0

        self.specfilehandle = None
        self.labels = None
        self.datablock = None

    @classmethod
    def initialize_esrf_ebs(cls):
        tmp =  at_data()

        tmp.at_file = '/scisoft/users/srio/Working/ATCOLLAB/esrf2_wiggler.spec'
        tmp.epsilon_x = 134e-12
        tmp.epsilon_y = 5e-12
        tmp.delta = 0.001,

        tmp.specfilehandle = specfile.Specfile(tmp.at_file)
        tmp.labels = numpy.array(tmp.specfilehandle[0].alllabels().copy())
        tmp.datablock = tmp.specfilehandle[0].data().copy()
        return tmp

    def get_data_with_label(self,label=""):

        igood = numpy.array(numpy.where(self.labels == label))
        if igood.size > 0:
            tmp = self.datablock[igood[0],:].copy()
            tmp.shape = -1
            return tmp

    def get_moments_at_sw(self,sw_value=0.0):
        sw = self.get_data_with_label("s-sW[m]")
        out = []
        for label in ["<xx>","<xxp>","<xpxp>","<yy>","<yyp>","<ypyp>"]:
            f1 = self.get_data_with_label(label)
            tmp =  numpy.interp(sw_value,sw,f1)
            out.append(tmp)
        return numpy.array(out)


    def get_moments_from_twiss(self,twiss=[0,0,0,0,0,0,0,0],dispersion=[0.0,0.0,0.0,0.0],epsilonx=134e-12,epsilony=5e-12,delta=0.001):
        mm = numpy.zeros(6)
        alphax = twiss[0]
        betax = twiss[1]
        gammax = twiss[2]

        alphay = twiss[3]
        betay = twiss[4]
        gammay = twiss[5]

        etax  = dispersion[0]
        etapx = dispersion[1]
        etay  = dispersion[2]
        etapy = dispersion[3]

        xx  =   betax  * epsilonx + (etax * delta)**2
        xxp = - alphax * epsilonx +  etax * etapx * delta**2
        xpxp =  gammax * epsilonx + (etapx * delta)**2

        yy  =   betay  * epsilony + (etay * delta)**2
        yyp = - alphay * epsilony +  etay * etapy * delta**2
        ypyp =  gammay * epsilony + (etapy * delta)**2
        return numpy.array([xx,xxp,xpxp,yy,yyp,ypyp])


    def get_twiss_at_sw(self,sw_value=0.0,propagate=0.0):
        sw = self.get_data_with_label("s-sW[m]")
        out = []
        for label in ["alphaX","betaX","gammaX","alphaY","betaY","gammaY"]:
            f1 = self.get_data_with_label(label)
            tmp =  numpy.interp(sw_value,sw,f1)
            out.append(tmp)
        out = numpy.array(out)

        if propagate == 0.0:
            return out
        else:
            alphax = out[0]
            betax = out[1]
            gammax = out[2]
            alphay = out[3]
            betay = out[4]
            gammay = out[5]
            out2 = numpy.zeros(6)
            out2[0] = alphax - gammax * propagate
            out2[1] = betax - 2 * alphax * propagate + gammax * propagate**2
            out2[2] = gammax
            out2[3] = alphay - gammay * propagate
            out2[4] = betay - 2 * alphay * propagate + gammay * propagate**2
            out2[5] = gammay
            return out2

    def get_dispersion_at_sw(self,sw_value=0.0,propagated=0.0):
        sw = self.get_data_with_label("s-sW[m]")
        out = []
        for label in ["etaX","etaXp","etaY","etaYp"]:
            f1 = self.get_data_with_label(label)
            tmp =  numpy.interp(sw_value,sw,f1)
            out.append(tmp)
        out = numpy.array(out)
        if propagated == 0.0:
            return out
        else:
            etax = out[0]
            etaxp = out[1]
            etay = out[2]
            etayp = out[3]
            out2 = numpy.zeros(4)
            out2[0] = etax + propagated * etaxp
            out2[1] = etaxp
            out2[2] = etay + propagated * etayp
            out2[3] = etayp
            return out2

    def get_sigmas_at_sw(self,sw_value=0.0):
        return 1e6*numpy.sqrt(numpy.abs(self.get_moments_at_sw(sw_value=sw_value)))

    def get_h_waists(self):
        #returns an array with the indices (float) where the zeroes are found (linear interpolation used)
        a = self.get_data_with_label("alphaX")
        return get_roots(a)

    def get_v_waists(self):
        #returns an array with the indices (float) where the zeroes are found (linear interpolation used)
        a = self.get_data_with_label("alphaY")
        return get_roots(a)

    def get_data_with_indices(self,indices=[1.1,2.2],label=""):
        x = self.get_data_with_label(label)
        i1 = (numpy.floor(indices)).astype(int)
        i2 =  (numpy.ceil(indices)).astype(int)
        # print(">><<",index,i1,i2)
        val = x[i1] + (indices - i1) * (x[i2] - x[i1])
        return val

def test_get_wiggler_center(position_center_wiggler=13.8379):
        ebs = at_data.initialize_esrf_ebs()
        cc = numpy.interp(0.0,ebs.get_data_with_label("s-sW[m]"),ebs.get_data_with_label("s[m]"))
        assert abs(position_center_wiggler - cc)<1e-3, 'Quick verification of wiggler center value'
        print("Absolute position of the wiggler center: %f m"%(cc))

def test_get_parameters_at_sw(sw_value=0.0):
    ebs = at_data.initialize_esrf_ebs()
    #
    # get values at given position
    #

    # print("======================================================================")
    # s_value = 13.8379
    # sw_value = s_value - position_center_wiggler  # -0.651 # 0.327 # 0.0 # -0.651
    #
    # print("\nCenter of wiggler at s=%f m"%(position_center_wiggler))
    print("\nInterpolated data at sw=%f m"%(sw_value))
    # print("\nAbsolute s=%f m"%(sw_value+position_center_wiggler))
    mm = ebs.get_moments_at_sw(sw_value)
    tt = ebs.get_twiss_at_sw(sw_value)
    ss = ebs.get_sigmas_at_sw(sw_value)

    print("<xx>:   %12g um ,   <xxp>: %12g   , <xpxp>: %12g"%(mm[0],mm[1],mm[2]))
    print("<zz>:   %12g um ,   <yyp>: %12g   , <ypyp>: %12g"%(mm[3],mm[4],mm[5]))
    print("sigmaX: %12.2f um , rhoX:%12g, sigmaX':%12.2f"%(ss[0],ss[1],ss[2]))
    print("sigmaY: %12.2f um , rhoY:%12g, sigmaY':%12.2f"%(ss[3],ss[4],ss[5]))
    print("alphaX: %12g      , betaX:  %12g,  gammaX: %12g"%(tt[0],tt[1],tt[2]))
    print("alphaY: %12g      , betaY:  %12g,  gammaY: %12g"%(tt[3],tt[4],tt[5]))

    assert abs(tt[0] - (-2.01801))<1e-2, 'Quick verification of alphaX'
    assert abs(tt[1] - (1.81379))<1e-2, 'Quick verification of betaX'
    assert abs(ss[3] - (3.58))<1e-2, 'Quick verification of sigmaY'
    assert abs(ss[5] - (3.04))<1e-2, 'Quick verification of sigmaYp'

def latex_table_twiss(file_out="table_twiss.txt"):
    #
    # create latex table with twiss parameters for different magnets
    #

    ebs = at_data.initialize_esrf_ebs()
    f = open(file_out,'w')
    devices = ["DQ2C\_2","short ID","QF8D","DQ1D"]
    s_values = [13.3872,13.8379,14.4302,15.0342]


    for i,s_value in enumerate(s_values):
        sw_value = s_value - position_center_wiggler
        tt = ebs.get_twiss_at_sw(sw_value)
        dd = ebs.get_dispersion_at_sw(sw_value)
        f.write("%15s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %d %s \n"%(
            devices[i],s_value,tt[0],tt[3],tt[1],tt[4],tt[2],tt[5],dd[0],dd[2],r"\\"))
    f.close()
    print("File written to disk: %s"%file_out)

def latex_table_waists(file_out="table_waists.txt"):
    #
    # create latex table with twiss parameters for waists
    #

    ebs = at_data.initialize_esrf_ebs()


    #
    # horizontal waist
    #
    h_idx = ebs.get_h_waists()
    h_sW = ebs.get_data_with_indices(h_idx,"s-sW[m]")
    tmp = get_roots(h_sW)
    wH_upstream = h_sW[int(tmp)]
    wH_downstream = h_sW[1+int(tmp)]

    #
    # vertical waist
    #
    idx = ebs.get_v_waists()
    sW = ebs.get_data_with_indices(idx,"s-sW[m]")
    tmp = get_roots(sW)
    wV_upstream = sW[int(tmp)]
    wV_downstream = sW[1+int(tmp)]

    #
    # create latex table with twiss parameters for waists
    #


    f = open(file_out,'w')

    tt = ebs.get_twiss_at_sw(wH_upstream)
    ss = ebs.get_sigmas_at_sw(wH_upstream)
    f.write("%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f &  %4.1f %s \n"%(
          r"H upstream  ($\beta$ min)",wH_upstream,wH_upstream+position_center_wiggler,tt[0],tt[1],ss[0],ss[2],ss[0]*ss[2],r"\\"))

    tt = ebs.get_twiss_at_sw(0.0)
    ss = ebs.get_sigmas_at_sw(0.0)
    f.write("%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %4.1f %s \n"%(
          r"H (center)",0,0+position_center_wiggler,tt[0],tt[1],ss[0],ss[2],ss[0]*ss[2],r"\\"))

    tt = ebs.get_twiss_at_sw(wH_downstream)
    ss = ebs.get_sigmas_at_sw(wH_downstream)
    f.write("%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %4.1f %s \n"%(
          r"H downstream  ($\beta$ max)",wH_downstream,wH_downstream+position_center_wiggler,tt[0],tt[1],ss[0],ss[2],ss[0]*ss[2],r"\\"))




    tt = ebs.get_twiss_at_sw(wV_upstream)
    ss = ebs.get_sigmas_at_sw(wV_upstream)
    f.write("%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %4.1f %s \n"%(
          r"V upstream  ($\beta$ max)",wV_upstream,wV_upstream+position_center_wiggler,tt[3],tt[4],ss[3],ss[5],ss[3]*ss[5],r"\\"))

    tt = ebs.get_twiss_at_sw(0.0)
    ss = ebs.get_sigmas_at_sw(0.0)
    f.write("%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %4.1f %s \n"%(
          r"V (center)",0,0+position_center_wiggler,tt[3],tt[4],ss[3],ss[5],ss[3]*ss[5],r"\\"))

    tt = ebs.get_twiss_at_sw(wV_downstream)
    ss = ebs.get_sigmas_at_sw(wV_downstream)
    f.write("%s & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %4.1f %s \n"%(
          r"V downstream  ($\beta$ min)",wV_downstream,wV_downstream+position_center_wiggler,tt[3],tt[4],ss[3],ss[5],ss[3]*ss[5],r"\\"))

    print("File written to disk: %s"%file_out)

def find_fake_waists(do_plot=0,eps_file=""):
    ebs = at_data.initialize_esrf_ebs()

    # create arrays extended the interval where the two closer real waists are found
    xx = numpy.linspace(-0.651*1.5,0.286*4,100)

    sigmax1 = numpy.zeros(xx.size)
    sigmax2 = numpy.zeros(xx.size)
    sigmay1 = numpy.zeros(xx.size)
    sigmay2 = numpy.zeros(xx.size)

    sigmaxp1 = numpy.zeros(xx.size)
    sigmaxp2 = numpy.zeros(xx.size)
    sigmayp1 = numpy.zeros(xx.size)
    sigmayp2 = numpy.zeros(xx.size)

    rhox1 = numpy.zeros(xx.size)
    rhox2 = numpy.zeros(xx.size)
    rhoy1 = numpy.zeros(xx.size)
    rhoy2 = numpy.zeros(xx.size)

    # fill arrays (real and propagated in absence of magnetic field)
    for i,x in enumerate(xx):
        mm = ebs.get_moments_at_sw(x)
        tt = ebs.get_twiss_at_sw(0.0,propagate=x)
        dd = ebs.get_dispersion_at_sw(0.0,propagated=x)
        mm1 = ebs.get_moments_from_twiss(tt,dd)
        sigmax1[i] = mm[0]
        sigmax2[i] = mm1[0]
        sigmay1[i] = mm[3]
        sigmay2[i] = mm1[3]

        rhox1[i] = mm[1]
        rhox2[i] = mm1[1]
        rhoy1[i] = mm[4]
        rhoy2[i] = mm1[4]

        sigmaxp1[i] = mm[2]
        sigmaxp2[i] = mm1[2]
        sigmayp1[i] = mm[5]
        sigmayp2[i] = mm1[5]

    # find the fake waists
    rrH = get_roots(rhox2)
    rrV = get_roots(rhoy2)
    rrH = rrH[0]
    rrV = rrV[0]
    fakeHwaist = interpolate_at_indices(xx,rrH)
    fakeVwaist = interpolate_at_indices(xx,rrV)

    print("H waists (idx, position, exact position, value, value at exact position): ",
          rrH,xx[int(rrH)],fakeHwaist, rhox2[int(rrH)],interpolate_at_indices(rhox2,rrH))
    print("V waists (idx, position, exact position, value, value at exact position): ",
          rrV,xx[int(rrV)],fakeVwaist, rhoy2[int(rrV)],interpolate_at_indices(rhoy2,rrV))

    print("H fake waist position: ",fakeHwaist," rho: ",interpolate_at_indices(rhox2,rrH))
    print("V fake waist position: ",fakeVwaist," rho: ",interpolate_at_indices(rhoy2,rrV))



    fake_sigma_h =  1e6*numpy.sqrt(interpolate_at_indices(sigmax2,rrH))
    fake_sigmap_h = 1e6*numpy.sqrt(interpolate_at_indices(sigmaxp2,rrH))
    fake_rho_h =    1e6*numpy.sqrt(numpy.abs(interpolate_at_indices(rhox2,rrH)))
    print("sigma, rho, sigmap, sigma*sigmap at H at waist: ",
        fake_sigma_h,fake_rho_h,fake_sigmap_h,fake_sigma_h*fake_sigmap_h)

    fake_sigma_v =  1e6*numpy.sqrt(interpolate_at_indices(sigmay2,rrV))
    fake_sigmap_v = 1e6*numpy.sqrt(interpolate_at_indices(sigmayp2,rrV))
    fake_rho_v =    1e6*numpy.sqrt(numpy.abs(interpolate_at_indices(rhoy2,rrV)))
    print("sigma, rho, sigmap, sigma*sigmap at V at waist: ",
        fake_sigma_v,fake_rho_v,fake_sigmap_v,fake_sigma_v*fake_sigmap_v)

    # make plot
    if do_plot:
        fig = plt.figure(1)
        axes = fig.add_subplot(211)
        # plt.plot(ebs.get_data_with_label("s-sW[m]"),ebs.get_data_with_label(label))
        plt.plot(xx,rhox1,label="real")
        plt.plot(xx,rhox2,label="propagated")
        plt.title("")
        plt.xlabel("$s_w$ [m]")
        plt.ylabel("$<xx'>^{1/2}$")
        plt.legend()
        # axes.set_xlim([-1,1])
        plt.grid()

        axes = fig.add_subplot(212)
        # plt.plot(ebs.get_data_with_label("s-sW[m]"),ebs.get_data_with_label(label))
        plt.plot(xx,rhoy1,label="real")
        plt.plot(xx,rhoy2,label="propagated")
        plt.title("")
        plt.xlabel("$s_w$ [m]")
        plt.ylabel("$<yy'>^{1/2}$")
        plt.legend()
        # axes.set_xlim([-1,1])
        plt.grid()
        if eps_file != "":
            plt.savefig(eps_file, format='eps', dpi=1000)
        plt.show()

    return fakeHwaist,fakeVwaist

def test_waists():

    ebs = at_data.initialize_esrf_ebs()
    #
    # horizontal waist
    #
    idx = ebs.get_h_waists()
    sW = ebs.get_data_with_indices(idx,"s-sW[m]")
    tmp = get_roots(sW)
    wH_upstream = sW[int(tmp)]
    wH_downstream = sW[1+int(tmp)]

    s =  ebs.get_data_with_indices(idx,"s[m]")
    t1 = ebs.get_data_with_indices(idx,"alphaX")
    t2 = ebs.get_data_with_indices(idx,"betaX")
    t3 = ebs.get_data_with_indices(idx,"gammaX")
    m1 = ebs.get_data_with_indices(idx,"SigX[um]")
    m2 = ebs.get_data_with_indices(idx,"<xxp>")
    m3 = ebs.get_data_with_indices(idx,"SigXp[urad]")
    #
    print("========== Horizontal waists ==========")
    for i in range(idx.size):
        print("i:%3d, idx:%9.3f, s:%6.3f, sW:%6.3f, a:%6.3f, b:%6.3f, g:%5.2f, sigmaX:%g, rhoX:%5.4f, sigmaY:%5.2f"%(
            i,idx[i],s[i],sW[i],t1[i],t2[i],t3[i],
            numpy.sqrt(m1[i]),
            1e6*numpy.sqrt(numpy.abs(m2[i])),
            numpy.sqrt(m3[i])) )


    #
    # vertical waist
    #
    idx = ebs.get_v_waists()
    sW = ebs.get_data_with_indices(idx,"s-sW[m]")
    tmp = get_roots(sW)
    wV_upstream = sW[int(tmp)]
    wV_downstream = sW[1+int(tmp)]


    s =  ebs.get_data_with_indices(idx,"s[m]")
    t1 = ebs.get_data_with_indices(idx,"alphaY")
    t2 = ebs.get_data_with_indices(idx,"betaY")
    t3 = ebs.get_data_with_indices(idx,"gammaY")
    m1 = ebs.get_data_with_indices(idx,"SigY[um]")
    m2 = ebs.get_data_with_indices(idx,"<yyp>")
    m3 = ebs.get_data_with_indices(idx,"SigYp[urad]")
    print(m3)
    print("========== Vertical waists ==========")
    for i in range(idx.size):
        print("i:%3d, idx:%9.3f, s:%6.3f, sW:%6.3f, a:%6.3f, b:%6.3f, g:%5.2f, sigmaY:%g, rhoY:%5.4f, sigmaY:%5.2f"%(
            i,idx[i],s[i],sW[i],t1[i],t2[i],t3[i],
            numpy.sqrt(m1[i]),
            1e6*numpy.sqrt(numpy.abs(m2[i])),
            numpy.sqrt(m3[i])) )


    assert abs(wH_upstream - (-0.651))<1e-3, 'Quick verification of position for H upstream waist '
    assert abs(wV_downstream - 0.286)<1e-3, 'Quick verification of position for V downstream waist '

def test_moments():
    #
    # values at wiggler center
    #
    ebs = at_data.initialize_esrf_ebs()

    sw_value = 0.0
    mm = ebs.get_moments_at_sw(sw_value)

    tt = ebs.get_twiss_at_sw(sw_value)
    dd = ebs.get_dispersion_at_sw(sw_value)
    mmtt = ebs.get_moments_from_twiss(tt,dd)

    print("Moments:                           ",mm)
    print("Moments from Twiss and Dispersion: ",mmtt)
    assert abs(numpy.abs(mm).sum() - numpy.abs(mmtt).sum())<1e-3, 'Quick verification moments consistency '

if __name__ == "__main__":

    position_center_wiggler = 13.8379
    #
    # some tests
    #
    test_waists()
    test_moments()
    test_get_wiggler_center(position_center_wiggler=position_center_wiggler)
    test_get_parameters_at_sw(sw_value=0.0)

    #
    # create latex tables for report
    #
    #latex_table_twiss(file_out="/users/srio/Working/rt/ESRF-new-lattice/REPORT-SHORT-WIGGLERS/table_twiss.txt")
    #latex_table_waists(file_out="/users/srio/Working/rt/ESRF-new-lattice/REPORT-SHORT-WIGGLERS/table_waists.txt")
    latex_table_twiss(file_out="table_twiss.txt")
    latex_table_waists(file_out="table_waists.txt")


    #
    #  find fake waists
    #
    #eps_file = "/users/srio/Working/rt/ESRF-new-lattice/REPORT-SHORT-WIGGLERS/GRAPHICS/fakewaists.eps"
    eps_file = "fakewaists.eps"
    (fakeHwaist,fakeVwaist) = find_fake_waists(do_plot=1,eps_file=eps_file)
