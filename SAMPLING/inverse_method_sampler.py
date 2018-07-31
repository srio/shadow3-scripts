

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
        return s1.get_sampled(cdf_rand_array)

    def get_n_sampled_points_and_histogram(self,npoints,bins=51,range=None):
        cdf_rand_array = numpy.random.random(npoints)
        return s1.get_sampled_and_histogram(cdf_rand_array,bins=bins,range=range)

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

if __name__ == "__main__":

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
