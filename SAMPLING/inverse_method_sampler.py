

import numpy

class Sampler1D(object):

    def __init__(self,pdf,pdf_x=None):

        self._pdf = pdf
        if pdf_x is None:
            self._set_default_pdf_x()
        else:
            self._pdf_x = pdf_x
        self._cdf = self._cdf_calculate()

    def pdf(self):
        return self._pdf

    def cdf(self):
        return self._cdf

    def abscissas(self):
        return self._pdf_x

    def get_sampled(self,random_in_0_1):
        y = numpy.array(random_in_0_1)

        if y.size > 1:
            x_rand_array = numpy.zeros_like(random_in_0_1)
            for i,cdf_rand in enumerate(random_in_0_1):
                ival,idelta,pendent = self._get_index(cdf_rand)
                x_rand_array[i] = self._pdf_x[ival] + idelta
            return x_rand_array
        else:
            ival,idelta,pendent = self._get_index(random_in_0_1)
            return self._pdf_x[ival] + idelta

    def get_sampled_and_histogram(self,random_in_0_1,bins=51,range=None):
        s1 = self.get_sampled(random_in_0_1)
        if range is None:
            range = [self._pdf_x.min(),self._pdf_x.max()]
        h,h_edges = numpy.array(numpy.histogram(s1,bins=bins,range=range))
        h_center = h_edges[0:-1] + 0.5*numpy.diff(h_edges)
        return s1,h,h_center

    def get_n_sampled_points(self,npoints):
        cdf_rand_array = numpy.random.random(npoints)
        return self.get_sampled(cdf_rand_array)

    def get_n_sampled_points_and_histogram(self,npoints,bins=51,range=None):
        cdf_rand_array = numpy.random.random(npoints)
        return self.get_sampled_and_histogram(cdf_rand_array,bins=bins,range=range)

    def _set_default_pdf_x(self):
        self._pdf_x = numpy.arange(self._pdf.size)

    def _cdf_calculate(self):
        cdf = numpy.cumsum(self._pdf)
        cdf -= cdf[0]
        cdf /= cdf.max()
        return cdf

    def _get_index(self,edge):
        ix = numpy.where(self._cdf > edge)[0][0]
        if ix > 0:
            ix -= 1
        if ix == (self._cdf.size-1):
            pendent = 0.0
            delta = 0.0
        else:
            pendent = self._cdf[ix+1] - self._cdf[ix]
            delta = (edge - self._cdf[ix]) / pendent
        return ix,delta,pendent

def test_1d():

    from srxraylib.plot.gol import plot

    x0=0.0
    sigma=2.0
    x = numpy.linspace(-10,10,51)
    y = numpy.exp(- (x-x0)**2 / 2 / sigma**2)

    y[0:21] = 100.0
    y[21:31] = 4.0
    y[31:41] = 5.0
    y[41:51] = 10.0


    s1 = Sampler1D(y,x)

    plot(s1.abscissas(),s1.pdf(),title="pdf")
    plot(s1.abscissas(),s1.cdf(),title="cdf")


    # defingn random points
    cdf_rand_array = numpy.random.random(100000)
    sampled_points,h,hx = s1.get_sampled_and_histogram(cdf_rand_array)

    plot(numpy.arange(cdf_rand_array.size),sampled_points,title="sampled points")
    plot(hx,h/h.max(),s1.abscissas(),s1.pdf()/s1.pdf().max(),title="histogram",legend=["histo","data"])


    # defining N
    sampled_points,h,hx = s1.get_n_sampled_points_and_histogram(120000)
    plot(numpy.arange(120000),sampled_points,title="120000 sapled points")
    plot(hx,h/h.max(),s1.abscissas(),s1.pdf()/s1.pdf().max(),title="histogram",legend=["histo","data"])

#
#
#
#


class Sampler2D(object):

    def __init__(self,pdf,pdf_x0=None,pdf_x1=None):

        self._pdf = pdf
        if pdf_x0 is None:
            self._set_default_pdf_x0()
        else:
            self._pdf_x0 = pdf_x0

        if pdf_x1 is None:
            self._set_default_pdf_x1()
        else:
            self._pdf_x1 = pdf_x1

        self._cdf2,self._cdf1 = self._cdf_calculate()

    def pdf(self):
        return self._pdf

    def cdf(self):
        return self._cdf2,self._cdf1

    def abscissas(self):
        return self._pdf_x0,self._pdf_x1

    def get_sampled(self,random0,random1):
        y0 = numpy.array(random0)
        y1 = numpy.array(random1)

        if y0.size > 1:
            x0_rand_array = numpy.zeros_like(y0)
            x1_rand_array = numpy.zeros_like(y1)

            for i,cdf_rand0 in enumerate(y0):
                ival,idelta,pendent = self._get_index0(cdf_rand0)
                x0_rand_array[i] = self._pdf_x0[ival] + idelta


                ival1,idelta1,pendent1 = self._get_index1(y1[i],ival)
                x1_rand_array[i] = self._pdf_x1[ival1] + idelta1
            return x0_rand_array,x1_rand_array
        else:
            pass
            # ival,idelta,pendent = self._get_index(random_in_0_1)
            # return self._pdf_x[ival] + idelta

    # def get_sampled_and_histogram(self,random_in_0_1,bins=51,range=None):
    #     pass
    #     # s1 = self.get_sampled(random_in_0_1)
    #     # if range is None:
    #     #     range = [self._pdf_x.min(),self._pdf_x.max()]
    #     # h,h_edges = numpy.array(numpy.histogram(s1,bins=bins,range=range))
    #     # h_center = h_edges[0:-1] + 0.5*numpy.diff(h_edges)
    #     # return s1,h,h_center

    def get_n_sampled_points(self,npoints):
        cdf_rand_array0 = numpy.random.random(npoints)
        cdf_rand_array1 = numpy.random.random(npoints)
        return self.get_sampled(cdf_rand_array0,cdf_rand_array1)

    # def get_n_sampled_points_and_histogram(self,npoints,bins=51,range=None):
    #     pass
    #     # cdf_rand_array = numpy.random.random(npoints)
    #     # return s1.get_sampled_and_histogram(cdf_rand_array,bins=bins,range=range)

    def _set_default_pdf_x(self):
        self._pdf_x = numpy.arange(self._pdf.size[0])

    def _set_default_pdf_y(self):
        self._pdf_y = numpy.arange(self._pdf.size[1])

    def _cdf_calculate(self):

        # cdf = numpy.cumsum(self._pdf)
        # cdf -= cdf[0]
        # cdf /= cdf.max()
        # return cdf

        cdf2 = numpy.zeros_like(self._pdf)
        cdf1 = numpy.zeros(cdf2.shape[0])
        for i in range(cdf2.shape[0]):
            cdf2[i,:] = numpy.cumsum(self._pdf[i,:])
            cdf1[i] = cdf2[i,-1]

            cdf2[i,:] -= cdf2[i,:][0]
            cdf2[i,:] /= cdf2[i,:].max()

        cdf1 = numpy.cumsum(cdf1)
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
            delta = (edge - self._cdf1[ix]) / pendent
        return ix,delta,pendent

    def _get_index1(self,edge,index0):
        ix = numpy.where(self._cdf2[index0,:] > edge)[0][0]
        if ix > 0:
            ix -= 1
        if ix == (self._cdf2[index0,:].size-1):
            pendent = 0.0
            delta = 0.0
        else:
            pendent = self._cdf2[index0,(ix+1)] - self._cdf2[index0,ix]
            delta = (edge - self._cdf2[index0,ix]) / pendent
        return ix,delta,pendent


def test_2d():
    from scipy.ndimage import imread
    from srxraylib.plot.gol import plot, plot_image, plot_scatter

    image_data = imread("/Users/srio/OASYS1.1/shadow3-scripts/SAMPLING/test1.jpg",flatten=True)
    image_data = numpy.flip(image_data.T,1)
    print(image_data.min(),image_data.max())
    image_data = image_data.max() - image_data
    # plot_image(image_data,cmap='binary')

    x0 = numpy.arange(image_data.shape[0])
    x1 = numpy.arange(image_data.shape[1])

    print(image_data.shape)

    s2d = Sampler2D(image_data,x0,x1)

    # plot_image(s2d.pdf(),cmap='binary',title="pdf")

    cdf2,cdf1 = s2d.cdf()
    plot_image(cdf2,cmap='binary',title="cdf")
    # plot(s2d.abscissas()[0],s2d.cdf()[0][:,-1])
    plot(s2d.abscissas()[0],cdf1)

    x0s,x1s = s2d.get_n_sampled_points(100000)
    plot_scatter(x0s,x1s)



if __name__ == "__main__":

    # test_1d()

    test_2d()


