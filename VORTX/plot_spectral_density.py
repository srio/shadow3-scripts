import numpy
import h5py
from srxraylib.plot.gol import plot_image,plot_surface, plot
import matplotlib.pylab as plt


from silx.math.fit.functions import sum_gauss
from silx.math.fit import fittheories
from silx.math.fit.fitmanager import FitManager

def get_fwhm(x,h):
    # h = ticket['histogram_h']
    tt = numpy.where(h>=max(h)*0.5)
    if h[tt].size > 1:
        binSize = x[1]-x[0]
        fwhm = binSize*(tt[0][-1]-tt[0][0])
        # ticket['fwhm_coordinates_h'] = (ticket['bin_h_center'][tt[0][0]],ticket['bin_h_center'][tt[0][-1]])
    else:
        fwhm = None

    return fwhm

def fit_gaussian(x,z1,do_plot=True):

    fwhm = get_fwhm(x,z1)
    print(" Graph FWHM is: ",fwhm)
    p = [z1.max(),0,fwhm]

    fit = FitManager()
    fit.setdata(x=x, y=z1)
    fit.loadtheories(fittheories)
    fit.settheory('Gaussians')
    fit.estimate()
    fit.runfit()

    print("Searched parameters = %s" % p)
    print("Obtained parameters : ")
    dummy_list = []
    for param in fit.fit_results:
        print(param['name'], ' = ', param['fitresult'])
        dummy_list.append(param['fitresult'])
    print("chisq = ", fit.chisq)
    fwhm_txt = "FWHM of fit = %5.3f um"%(fit.fit_results[2]['fitresult'])

    z11 = sum_gauss(x, *dummy_list)

    if do_plot:
        plot(x,z1,x,z11,legend=["data","fit"],ylog=False)

    FWHM = fit.fit_results[2]['fitresult']

    return FWHM

def spectral_density(h5file_root="vx_id16a_A",up_to_mode=1099,do_plot=False):


    h5file = h5file_root+".h5"

    f = h5py.File(h5file,'r')

    arr1 = f["uptomode%04d/SpectralDensity/image_data"%up_to_mode].value.T
    x    = f["uptomode%04d/SpectralDensity/axis_x"%up_to_mode].value
    y    = f["uptomode%04d/SpectralDensity/axis_y"%up_to_mode].value
    f.close()

    # print(arr1[0,0])
    fig = plot_image(arr1,1e3*x,1e3*y,cmap='jet',show=0) # ,xrange=[-50,50],yrange=[-50,50])

    # fig = plt.figure()
    # plt.imshow(arr1.T,cmap='jet',extent=[1e3*x[0],1e3*x[-1],1e3*y[0],1e3*y[-1]],interpolation='none',aspect='auto')
    # #


    # if up_to_mode == 0:
    # plt.colorbar(orientation="horizontal", pad=0.2)
    # else:
    #     ax = fig.gca()
    #     ax.set_xlabel("X [$\mu$m]")
    #     ax.set_ylabel("Y [$\mu$m]")

    plt.savefig("/tmp/spectral_density.png",dpi=600)

    plt.show()

    z1 = arr1.sum(axis=1)
    fwhm_x = fit_gaussian(x,z1,do_plot=do_plot)

    z0 = arr1.sum(axis=0)
    fwhm_y = fit_gaussian(y,z0,do_plot=do_plot)

    print("FWHMs [um]: ",1e3*fwhm_x,1e3*fwhm_y)




if __name__ == "__main__":

    spectral_density(h5file_root="vx_id16a_A",up_to_mode=1099,do_plot=True)
