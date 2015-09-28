#
# overwrite the energy in a SHADOW source with values sampled from 
# a spectrum (defined as an array)
#
# srio@esrf.eu  2015-07-15
#


import numpy
import Shadow


def sample_from_spectrum(spectrum,npoints=1000,plot=0):
    # normalize distribution function 
    spectrum[:,1] /= spectrum[:,1].sum()
    #calculate the cumulative distribution function
    a_cdf = numpy.cumsum(spectrum[:,1])
    # create randomly distributed npoints 
    rd = numpy.random.rand(npoints)
    # evaluate rd following the inverse of cdf
    sampled_energies = numpy.interp(rd,a_cdf,spectrum[:,0])
    #plots
    if plot:
        import matplotlib.pylab as plt
        plt.subplot(111)
        plt.figure(1)
        plt.xlabel("photon energy /eV")
        plt.ylabel("normalized spectrum")
        plt.plot(a[:,0],a[:,1])
        
        plt.figure(2)
        plt.xlabel("a_cdf")
        plt.ylabel("photon energy / eV")
        plt.plot(a_cdf,a[:,0])
        
        plt.figure(3)
        plt.xlabel("index")
        plt.ylabel("photon energy / eV")
        plt.plot(numpy.linspace(0,npoints-1,npoints), sampled_energies)

        plt.show()
    
    return sampled_energies


if __name__ == "__main__":
    
    #define distribution function (spectrum)
    a = numpy.array([ \
        [1000.00,0.00000 ],\
        [2000.00,0.00000 ],\
        [3000.00,0.00000 ],\
        [4000.00,0.00000 ],\
        [5000.00,0.00000 ],\
        [6000.00,0.00000 ],\
        [7000.00,0.00000 ],\
        [8000.00,0.00000 ],\
        [9000.00,0.00000 ],\
        [10000.0,6.65281 ],\
        [11000.0,26.1984 ],\
        [12000.0,86.3888 ],\
        [13000.0,285.715 ],\
        [14000.0,642.834 ],\
        [15000.0,1605.83 ],\
        [16000.0,2957.33 ],\
        [17000.0,6046.30 ],\
        [18000.0,9838.37 ],\
        [19000.0,12738.8 ],\
        [20000.0,16130.6 ],\
        [21000.0,19788.5 ],\
        [22000.0,23865.6 ],\
        [23000.0,27854.9 ],\
        [24000.0,32146.1 ],\
        [25000.0,35882.3 ],\
        [26000.0,39455.6 ],\
        [27000.0,42370.8 ],\
        [28000.0,45281.9 ],\
        [29000.0,47336.9 ],\
        [30000.0,49444.8 ],\
        [31000.0,50532.1 ],\
        [32000.0,51625.8 ],\
        [33000.0,52335.6 ],\
        [34000.0,53043.7 ],\
        [35000.0,53127.2 ],\
        [36000.0,53049.9 ],\
        [37000.0,52808.1 ],\
        [38000.0,52491.1 ],\
        [39000.0,52040.4 ],\
        [40000.0,51575.2 ],\
        [41000.0,50810.7 ],\
        [42000.0,50050.4 ],\
        [43000.0,49051.7 ],\
        [44000.0,48092.7 ],\
        [45000.0,47061.9 ],\
        [46000.0,45941.3 ],\
        [47000.0,45131.6 ],\
        [48000.0,44310.4 ],\
        [49000.0,42984.9 ],\
        [50000.0,41634.4 ],\
        [51000.0,40430.2 ],\
        [52000.0,39222.2 ],\
        [53000.0,38071.3 ],\
        [54000.0,36896.7 ],\
        [55000.0,35945.7 ],\
        [56000.0,34981.7 ],\
        [57000.0,42099.2 ],\
        [58000.0,49257.8 ],\
        [59000.0,55207.9 ],\
        [60000.0,61086.8 ],\
        [61000.0,46152.4 ],\
        [62000.0,31065.0 ],\
        [63000.0,28571.1 ],\
        [64000.0,26042.6 ],\
        [65000.0,25171.6 ],\
        [66000.0,24200.1 ],\
        [67000.0,28353.4 ],\
        [68000.0,32520.8 ],\
        [69000.0,27392.2 ],\
        [70000.0,22032.3 ],\
        [71000.0,19236.8 ],\
        [72000.0,16385.0 ],\
        [73000.0,15587.5 ],\
        [74000.0,14719.0 ],\
        [75000.0,14176.4 ],\
        [76000.0,13719.6 ],\
        [77000.0,13059.6 ],\
        [78000.0,12319.1 ],\
        [79000.0,11773.1 ],\
        [80000.0,10935.6 ],\
        [81000.0,10549.2 ],\
        [82000.0,9802.11 ],\
        [83000.0,9181.41 ],\
        [84000.0,8716.42 ],\
        [85000.0,8071.86 ],\
        [86000.0,7430.28 ],\
        [87000.0,6998.92 ],\
        [88000.0,6568.78 ],\
        [89000.0,5975.87 ],\
        [90000.0,5382.43 ],\
        [91000.0,4901.83 ],\
        [92000.0,4425.00 ],\
        [93000.0,3795.52 ],\
        [94000.0,3258.13 ],\
        [95000.0,2786.15 ],\
        [96000.0,2312.45 ],\
        [97000.0,1879.51 ],\
        [98000.0,1430.59 ],\
        [99000.0,920.485 ],\
        [100000.,414.051 ],\
        [101000.,95.4982 ],\
        [102000.,0.00000 ],\
        [103000.,0.00000 ] ] )
    
    print("Energy from %f eV to %f eV"%(a[0,0],a[-1,0]))
    #
    # create SHADOW source 
    #
    
    src = Shadow.Source()
    
    beam = Shadow.Beam()
    beam.genSource(src)
    
    npoints = src.NPOINT
    sampled_energies = sample_from_spectrum(a,npoints=npoints,plot=1)
    
    #calculate wavevectors = 2 pi / lambda
    codata_h = numpy.array(6.62606957e-34)
    codata_ec = numpy.array(1.602176565e-19)
    codata_c = numpy.array(299792458.0)
    A2EV = 2.0*numpy.pi/(codata_h*codata_c/codata_ec*1e2) #A2EV = 50676.89919
    sampled_wavevectors = sampled_energies * A2EV #2 pi / wavelength in cm^-1

    #overwrite the energy column with new values
    beam.rays[:,10] = sampled_wavevectors 
    
    
    # calculate the histogram for energy of the new source
    import Shadow.ShadowTools 
    Shadow.ShadowTools.histo1(beam,11,nbins=75)
    
