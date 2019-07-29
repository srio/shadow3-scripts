# Hi
# Manuel,
#
# Here
# 's a python version of the fine determination of the FWHM. It allows to be much more precise that using argmax, for instance I have these results:
# coarse
# fhwm: 122.000
# px
# actual
# fwhm: 123.456
# px
# fine
# fhwm: 123.457
# px
#
# It
# 's useful when running loops and simulation.
#
# In
# the
# case
# of
# a
# random
# gaussian
# process, the
# issue
# indeed is that
# Cramer - Rao
# lower
# bound
# comes
# for the kill

# awojdyla@ lbl.gov, July 2019
import numpy as np

def get_fwhm_coarse(signal):

    # half maximum (could be changed)
    thr = np.amax(signal / 2)

    bounds = np.argwhere(signal > thr)

    # here are the coarse estimate of threshold crossing (like in histo1)
    ixl_e = np.amin(bounds)
    ixr_e = np.amax(bounds)

    fwhm_coarse_px = np.abs(ixl_e - ixr_e)

    return ixl_e,ixr_e,fwhm_coarse_px

def get_fwhm_fine(signal):

    # half maximum (could be changed)
    thr = np.amax(signal / 2)

    bounds = np.argwhere(signal > thr)

    # here are the coarse estimate of threshold crossing (like in histo1)
    ixl_e = np.amin(bounds);
    ixr_e = np.amax(bounds);

    # refine the threasold crossing estimate using
    # explicit linear interpolation

    # left edge
    if ixl_e > 0:  # make sure there is a left edge
        xl = ixl_e - (signal[ixl_e] - thr) / (signal[ixl_e] - signal[ixl_e - 1]);
    else:  # otherwise, pupulate missing edge as NaNs
        xl = np.nan

    # right edge
    if ixr_e < np.size(signal):
        xr = ixr_e - (signal[ixr_e] - thr) / (signal[ixr_e + 1] - signal[ixr_e]);
    else:
        xr = np.nan

    fwhm_fine_px = np.abs(xr - xl)

    return xl,xr,fwhm_fine_px


if __name__ == "__main__":
    # center of a gaussian
    mean_px = 0
    # fwhm of the gaussian
    fwhm_px = 123.456

    # number of points for the gaussian
    N_pts = 1000
    x_px = np.arange(-N_pts / 2, N_pts / 2)

    # define the gaussian
    sigma_x = fwhm_px / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-((x_px - mean_px) / (np.sqrt(2) * sigma_x)) ** 2) + (np.random.rand(N_pts) - 0.5) * 0.1

    # plot the gaussian
    from matplotlib import pyplot as plt

    plt.plot(x_px, g, "+")
    plt.plot(x_px, g)
    plt.xlabel('position [px]')
    plt.ylabel('intensity')


    # take the gaussian as an example for the prod
    # signal = g

    xl,xr,fwhm_coarse_px = get_fwhm_coarse(g)
    yl,yr = x_px[xl],x_px[xr]

    xl0,xr0,fwhm_fine_px = get_fwhm_fine(g)
    yl0 = np.interp(xl0,range(x_px.size),x_px)
    yr0 = np.interp(xr0, range(x_px.size), x_px)



    plt.plot((yl0, yr0), (0.5 * g.max(), 0.5 * g.max()))
    plt.plot((yl, yr), (0.5 * g.max(), 0.5 * g.max()))
    plt.show()
    print( 'coarse  fhwm: %1.3f px' % fwhm_coarse_px +
          '\nactual fwhm: %1.3f px' % fwhm_px +
          '\nfine   fhwm: %1.3f px' % fwhm_fine_px)

    print(xl,xr,xl0,xr0)