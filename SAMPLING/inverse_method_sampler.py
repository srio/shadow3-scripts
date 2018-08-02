import numpy

class Sampler1D(object):

    def __init__(self,pdf,pdf_x=None):

        self._pdf = pdf
        if pdf_x is None:
            self._set_default_x()
        else:
            self._x = pdf_x
        self._cdf = self._cdf_calculate()

    def pdf(self):
        return self._pdf

    def cdf(self):
        return self._cdf

    def abscissas(self):
        return self._x

    def get_sampled(self,random_in_0_1):
        y = numpy.array(random_in_0_1)

        if y.size > 1:
            x_rand_array = numpy.zeros_like(random_in_0_1)
            for i,cdf_rand in enumerate(random_in_0_1):
                ival,idelta,pendent = self._get_index(cdf_rand)
                x_rand_array[i] = self._x[ival] + idelta*(self._x[1]-self._x[0])
            return x_rand_array
        else:
            ival,idelta,pendent = self._get_index(random_in_0_1)
            return self._x[ival] + idelta*(self._x[1]-self._x[0])

    def get_n_sampled_points(self,npoints):
        cdf_rand_array = numpy.random.random(npoints)
        return self.get_sampled(cdf_rand_array)

    def get_sampled_and_histogram(self,random_in_0_1,bins=51,range=None):
        s1 = self.get_sampled(random_in_0_1)
        if range is None:
            range = [self._x.min(),self._x.max()]
        h,h_edges = numpy.array(numpy.histogram(s1,bins=bins,range=range))
        h_center = h_edges[0:-1] + 0.5*numpy.diff(h_edges)
        return s1,h,h_center

    def get_n_sampled_points_and_histogram(self,npoints,bins=51,range=None):
        cdf_rand_array = numpy.random.random(npoints)
        return self.get_sampled_and_histogram(cdf_rand_array,bins=bins,range=range)

    def _set_default_x(self):
        self._x = numpy.arange(self._pdf.size)

    def _cdf_calculate(self):
        cdf = numpy.cumsum(self._pdf)
        if (cdf[-1] - cdf[0]) > 1e-10:
            cdf -= cdf[0]
        cdf /= cdf[-1]
        return cdf

    def _get_index(self,edge):
        ix = numpy.where(self._cdf > edge)[0][0]
        if ix > 0:
            ix -= 1
        if ix == (self._cdf.size-1):
            pendent = 0.0
            delta = 0.0
        else:
            pendent = (self._cdf[ix+1] - self._cdf[ix])
            # if (ix+2 <= self._cdf.size-1): pendent = 0.5*(self._cdf[ix+2] - self._cdf[ix])
            delta = 0.0
            if numpy.abs(pendent) > 1e-10:
                delta = (edge - self._cdf[ix]) / pendent
        # delta *= 3./2
        # print(">>>",ix,ix+1,delta)
        # print("   ",self._x[ix],self._x[ix+1],self._x[ix]+delta*(self._x[ix+1]-self._x[ix]))
        # print("   ",self._cdf[ix],self._cdf[ix+1],self._cdf[ix]+delta*pendent,edge)
        return ix,delta,pendent

def test_1d_gaussian():

    from srxraylib.plot.gol import plot

    x0=0.0
    sigma=2.0
    x = numpy.linspace(-10,10,51)
    y = numpy.exp(- (x-x0)**2 / 2 / sigma**2)


    s1 = Sampler1D(y,x)

    # plot(s1.abscissas(),s1.pdf(),title="pdf")
    plot(s1.abscissas(),s1.cdf(),title="cdf")
    # for i in range(s1.abscissas().size):
    #     print(s1.abscissas()[i],s1.pdf()[i],s1.cdf()[i])


    # defingn random points
    cdf_rand_array = numpy.random.random(100000)
    sampled_points,h,hx = s1.get_sampled_and_histogram(cdf_rand_array)


    # plot(numpy.arange(cdf_rand_array.size),sampled_points,title="sampled points")
    plot(hx,h/h.max(),s1.abscissas(),s1.pdf()/s1.pdf().max(),title="histogram",legend=["histo","data"])


    # defining N
    # sampled_points,h,hx = s1.get_n_sampled_points_and_histogram(120000)
    # plot(numpy.arange(120000),sampled_points,title="120000 sapled points")
    # plot(hx,h/h.max(),s1.abscissas(),s1.pdf()/s1.pdf().max(),title="histogram",legend=["histo","data"])
    print("sampled points mean, stdev, step: ",sampled_points.mean(),sampled_points.std(),x[1]-x[0])
#
def test_1d_external_cdf():

    from srxraylib.plot.gol import plot

    x0=0.0
    sigma=2.0
    x = numpy.linspace(-10,10,51)
    y = numpy.exp(- (x-x0)**2 / 2 / sigma**2)


    s1 = Sampler1D(y,x)

    # overwrite _cdf
    from math import erf,sqrt

    cdf = numpy.zeros_like(s1.abscissas())
    for i in range(cdf.size):
        cdf[i] = (1.0 + erf(s1.abscissas()[i] / sigma / sqrt(2.0))) / 2.0
    s1._cdf = cdf


    plot(s1.abscissas(),s1.pdf(),title="pdf")
    plot(s1.abscissas(),s1.cdf(),title="cdf")
    for i in range(s1.abscissas().size):
        print(s1.abscissas()[i],s1.pdf()[i],s1.cdf()[i])


    # defingn random points
    cdf_rand_array = numpy.random.random(100000)
    sampled_points,h,hx = s1.get_sampled_and_histogram(cdf_rand_array)


    # plot(numpy.arange(cdf_rand_array.size),sampled_points,title="sampled points")
    plot(hx,h/h.max(),s1.abscissas(),s1.pdf()/s1.pdf().max(),title="histogram",legend=["histo","data"])


    # defining N
    # sampled_points,h,hx = s1.get_n_sampled_points_and_histogram(120000)
    # plot(numpy.arange(120000),sampled_points,title="120000 sapled points")
    # plot(hx,h/h.max(),s1.abscissas(),s1.pdf()/s1.pdf().max(),title="histogram",legend=["histo","data"])
    print("sampled points mean, stdev, step: ",sampled_points.mean(),sampled_points.std(),x[1]-x[0])
#

def test_1d_patological():

    from srxraylib.plot.gol import plot

    x0=0.0
    sigma=2.0
    x = numpy.linspace(-10,10,501)
    y = x*0.0

    # y[-1] = 100.0

    # print(">>>>>>>",x[200])
    # y[200] = 300.0


    y[0] = 4.0



    s1 = Sampler1D(y,x)

    plot(s1.abscissas(),s1.pdf(),title="pdf",marker="o")
    plot(s1.abscissas(),s1.cdf(),title="cdf",marker="o")


    # defingn random points
    cdf_rand_array = numpy.random.random(100000)
    sampled_points,h,hx = s1.get_sampled_and_histogram(cdf_rand_array)

    plot(numpy.arange(cdf_rand_array.size),sampled_points,title="sampled points")
    plot(hx,h/h.max(),s1.abscissas(),s1.pdf()/s1.pdf().max(),title="histogram",legend=["histo","data"])

def test_1d_single_point():


    x = numpy.linspace(-10,10,501)
    y = numpy.abs(x)

    s1 = Sampler1D(y,x)

    # from srxraylib.plot.gol import plot
    # plot(s1.abscissas(),s1.pdf(),title="pdf",marker="o")
    # plot(s1.abscissas(),s1.cdf(),title="cdf",marker="o")


    # defining random points
    cdf_rand_array = numpy.random.random(100000)
    sampled_points = s1.get_sampled(cdf_rand_array)

    sampled_point_200 = s1.get_sampled(cdf_rand_array[200])

    print(cdf_rand_array[200],sampled_points[200],sampled_point_200)
    assert(sampled_points[200] == sampled_point_200)



class Sampler2D(object):

    def __init__(self,pdf,x0=None,x1=None):

        self._pdf = pdf
        if x0 is None:
            self._set_default_x0()
        else:
            self._x0 = x0

        if x1 is None:
            self._set_default_x1()
        else:
            self._x1 = x1

        self._cdf2,self._cdf1 = self._cdf_calculate()

    def pdf(self):
        return self._pdf

    def cdf(self):
        return self._cdf2,self._cdf1

    def abscissas(self):
        return self._x0,self._x1

    def get_sampled(self,random0,random1):
        if random0.size != random1.size:
            raise Exception("Dimension of two random arrays must be equal.")

        if random0.size > 1:
            y0 = numpy.array(random0)
            y1 = numpy.array(random1)
        else:
            y0 = numpy.array([random0])
            y1 = numpy.array([random1])

        x0_rand_array = numpy.zeros_like(y0)
        x1_rand_array = numpy.zeros_like(y1)

        for i,cdf_rand0 in enumerate(y0):
            index0,idelta0,pendent0 = self._get_index0(cdf_rand0)
            x0_rand_array[i] = self._x0[index0] + idelta0*(self._x0[1]-self._x0[0])

            index1,idelta1,pendent1 = self._get_index1(y1[i],index0,delta0=idelta0)
            x1_rand_array[i] = self._x1[index1] + idelta1*(self._x1[1]-self._x1[0])

        if y0.size > 1:
            return x0_rand_array,x1_rand_array
        else:
            return x0_rand_array[0],x1_rand_array[0]

    def get_n_sampled_points(self,npoints):
        cdf_rand_array0 = numpy.random.random(npoints)
        cdf_rand_array1 = numpy.random.random(npoints)
        return self.get_sampled(cdf_rand_array0,cdf_rand_array1)

    # def get_n_sampled_points_and_histogram(self,npoints,bins=51,range=None):
    #     pass
    #     # cdf_rand_array = numpy.random.random(npoints)
    #     # return s1.get_sampled_and_histogram(cdf_rand_array,bins=bins,range=range)

    def _set_default_x(self):
        self._x = numpy.arange(self._pdf.size[0])

    def _set_default_pdf_y(self):
        self._pdf_y = numpy.arange(self._pdf.size[1])

    def _cdf_calculate(self):
        cdf2 = numpy.zeros_like(self._pdf)
        cdf1 = numpy.zeros(cdf2.shape[0])
        for i in range(cdf2.shape[0]):
            cdf2[i,:] = numpy.cumsum(self._pdf[i,:])
            cdf1[i] = cdf2[i,-1]

            if (cdf2[i,-1] - cdf2[i,0]) > 1e-10:
                cdf2[i,:] -= cdf2[i,:][0]
            cdf2[i,:] /= cdf2[i,:][-1]

        cdf1 = numpy.cumsum(cdf1)
        if (cdf1[-1] - cdf1[0]) > 1e-10:
            cdf1 -= cdf1[0]
        cdf1 /= cdf1.max()

        return cdf2,cdf1

    def _get_index0(self,edge):
        ix = numpy.where(self._cdf1 > edge)[0][0]
        if ix > 0:
            ix -= 1
        if ix == (self._cdf1.size-1):
            pendent = 0.0
            delta = 0.0
        else:
            pendent = self._cdf1[ix+1] - self._cdf1[ix]
            delta = 0.0
            if numpy.abs(pendent) > 1e-10:
                delta = (edge - self._cdf1[ix]) / pendent
        return ix,delta,pendent

    def _get_index1(self,edge,index0,delta0=None):
        ix = numpy.where(self._cdf2[index0,:] > edge)[0][0]
        if ix > 0:
            ix -= 1
        if ix == (self._cdf2[index0,:].size-1):
            pendent = 0.0
            delta = 0.0
        else:
            if delta0 is None:
                cdf2 = self._cdf2[index0,:]
            else:
                cdf2 = self._cdf2[index0,:] + delta0*(self._cdf2[index0+1,:]-self._cdf2[index0,:])
            pendent = cdf2[(ix+1)] - cdf2[ix]
            delta = 0.0
            if numpy.abs(pendent) > 1e-10:
                delta = (edge - cdf2[ix]) / pendent
        return ix,delta,pendent


def test_2d_image():
    from scipy.ndimage import imread
    from srxraylib.plot.gol import plot, plot_image, plot_scatter

    image_data = imread("test1.jpg",flatten=True)
    image_data = numpy.flip(image_data.T,1)
    image_data = image_data.max() - image_data


    x0 = numpy.arange(image_data.shape[0])
    x1 = numpy.arange(image_data.shape[1])

    s2d = Sampler2D(image_data,x0,x1)

    plot_image(s2d.pdf(),cmap='binary',title="pdf")

    cdf2,cdf1 = s2d.cdf()
    plot_image(cdf2,cmap='binary',title="cdf")
    plot(s2d.abscissas()[0],cdf1)

    x0s,x1s = s2d.get_n_sampled_points(100000)
    plot_scatter(x0s,x1s)

def test_2d_image_ascent():
    import scipy.misc
    from srxraylib.plot.gol import plot, plot_image, plot_scatter


    image_data = numpy.array(scipy.misc.ascent(),dtype=float)

    image_data = image_data.max() - numpy.rot90(numpy.rot90(numpy.rot90(image_data,1)))

    x0 = numpy.arange(image_data.shape[0])
    x1 = numpy.arange(image_data.shape[1])

    s2d = Sampler2D(image_data,x0,x1)

    plot_image(s2d.pdf(),cmap='binary',title="pdf")

    cdf2,cdf1 = s2d.cdf()
    plot_image(cdf2,cmap='binary',title="cdf")
    plot(s2d.abscissas()[0],cdf1)

    x0s,x1s = s2d.get_n_sampled_points(200000)
    plot_scatter(x0s,x1s)


def test_2d_gaussian():
    from srxraylib.plot.gol import plot, plot_image, plot_scatter

    x0 = numpy.linspace(-800.0,800,200)
    x1 = numpy.linspace(-350.0,350,100)

    X0 = numpy.outer(x0,numpy.ones_like(x1))
    X1 = numpy.outer(numpy.ones_like(x0),x1)

    sigmax0 = 150.4
    sigmax1 = 100.0
    angle_in_deg = 0.0

    angle = angle_in_deg * numpy.pi / 180.0

    a = (numpy.cos(angle))**2 / (2 * sigmax0**2) + (numpy.sin(angle))**2 / (2 * sigmax1**2)
    c = (numpy.sin(angle))**2 / (2 * sigmax0**2) + (numpy.cos(angle))**2 / (2 * sigmax1**2)
    b = -numpy.sin(2*angle) / (4 * sigmax0**2) + -numpy.sin(2*angle) / (4 * sigmax1**2)
    Z = numpy.exp( -1.0*(a*X0**2 + b*X0*X1 + c*X1**2))
    print(a,b,c,X0)

    # plot_image(Z,x0,x1)

    s2d = Sampler2D(Z,x0,x1)

    plot_image(s2d.pdf(),s2d.abscissas()[0],s2d.abscissas()[1],cmap='binary',title="pdf")

    cdf2,cdf1 = s2d.cdf()
    plot_image(cdf2,cmap='binary',title="cdf")
    # # plot(s2d.abscissas()[0],s2d.cdf()[0][:,-1])
    plot(s2d.abscissas()[0],cdf1)

    x0s,x1s = s2d.get_n_sampled_points(100000)
    plot_scatter(x0s,x1s,nbins=151)
    # plot_scatter(numpy.arange(x0s.size),x0s)




def test_2d_single_point():
    from scipy.ndimage import imread
    from srxraylib.plot.gol import plot, plot_image, plot_scatter

    x0 = numpy.linspace(-800.0,800,200)
    x1 = numpy.linspace(-350.0,350,100)

    X0 = numpy.outer(x0,numpy.ones_like(x1))
    X1 = numpy.outer(numpy.ones_like(x0),x1)

    sigmax0 = 150.4
    sigmax1 = 100.0
    angle_in_deg = 0.0

    angle = angle_in_deg * numpy.pi / 180.0

    a = (numpy.cos(angle))**2 / (2 * sigmax0**2) + (numpy.sin(angle))**2 / (2 * sigmax1**2)
    c = (numpy.sin(angle))**2 / (2 * sigmax0**2) + (numpy.cos(angle))**2 / (2 * sigmax1**2)
    b = -numpy.sin(2*angle) / (4 * sigmax0**2) + -numpy.sin(2*angle) / (4 * sigmax1**2)
    Z = numpy.exp( -1.0*(a*X0**2 + b*X0*X1 + c*X1**2))
    print(a,b,c,X0)

    # plot_image(Z,x0,x1)

    s2d = Sampler2D(Z,x0,x1)

    # plot_image(s2d.pdf(),s2d.abscissas()[0],s2d.abscissas()[1],cmap='binary',title="pdf")

    cdf2,cdf1 = s2d.cdf()
    # plot_image(cdf2,cmap='binary',title="cdf")
    # # plot(s2d.abscissas()[0],s2d.cdf()[0][:,-1])
    # plot(s2d.abscissas()[0],cdf1)


    rand_array0 = numpy.random.random(1000)
    rand_array1 = numpy.random.random(1000)

    x0s,x1s = s2d.get_sampled(rand_array0,rand_array1)
    x0s200,x1s200 = s2d.get_sampled(rand_array0[200],rand_array1[200])

    print(rand_array0[200],rand_array1[200],x0s[200],x1s[200],x0s200,x1s200)
    assert(x0s[200] == x0s200)
    assert(x1s[200] == x1s200)



if __name__ == "__main__":

    # test_1d_gaussian()

    test_1d_external_cdf()
    # test_1d_patological()

    # test_2d_image()

    # test_2d_image_ascent()

    # test_2d_gaussian()

    # test_1d_single_point()

    # test_2d_single_point()


    # from math import erf,sqrt
    # x = 10.0
    # sigma = 2.0
    # print((1.0 + erf(x / sigma / sqrt(2.0))) / 2.0)