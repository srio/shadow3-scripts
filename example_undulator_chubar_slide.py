"""

    This example shows the focusing of an ideal lens in 1:1 configuration
    for an undulator source sources (see main program at the bottom)


"""
import numpy

from srxraylib.waveoptics.wavefront2D  import Wavefront2D
from srxraylib.waveoptics.propagator2D import propagate_2D_fraunhofer
from srxraylib.waveoptics.propagator2D import propagate_2D_fresnel, propagate_2D_fresnel_convolution, propagate_2D_fresnel_srw

import scipy.constants as codata

# this is undulator block
from pySRU.ElectronBeam import ElectronBeam
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane
from pySRU.TrajectoryFactory import TrajectoryFactory, TRAJECTORY_METHOD_ANALYTIC,TRAJECTORY_METHOD_ODE
from pySRU.RadiationFactory import RadiationFactory,RADIATION_METHOD_NEAR_FIELD, RADIATION_METHOD_APPROX_FARFIELD
from pySRU.SourceUndulatorPlane import SourceUndulatorPlane

try:
    import srwlib
    SRWLIB_AVAILABLE = True
except:
    SRWLIB_AVAILABLE = False
    print("SRW is not available")


# global variables
codata_mee = numpy.array(codata.physical_constants["electron mass energy equivalent in MeV"][0])
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)

do_plot = True

if do_plot:
    from srxraylib.plot.gol import plot,plot_image,plot_table

#
# auxiliar functions
#

def line_fwhm(line):
    #
    #CALCULATE fwhm in number of abscissas bins (supposed on a regular grid)
    #
    tt = numpy.where(line>=max(line)*0.5)
    if line[tt].size > 1:
        # binSize = x[1]-x[0]
        FWHM = (tt[0][-1]-tt[0][0])
        return FWHM
    else:
        return -1


def propagation_to_image(wf,do_plot=do_plot,plot_title="Before lens",method='fft',
                            propagation_distance=30.0,defocus_factor=1.0,propagation_steps=1,show=1):


    method_label = "fresnel (%s)"%method
    print("\n#                                                             ")
    print("# near field fresnel (%s) diffraction and focusing  "%(method_label))
    print("#                                                             ")

    #                               \ |  /
    #   *                           | | |                      *
    #                               / | \
    #   <-------    d  ---------------><---------   d   ------->
    #   d is propagation_distance

    print("Incident intensity: ",wf.get_intensity().sum())

    # propagation downstream the lens to image plane
    for i in range(propagation_steps):
        if propagation_steps > 1:
            print(">>> Propagating step %d of %d; propagation_distance=%g m"%(i+1,propagation_steps,
                                                propagation_distance*defocus_factor/propagation_steps))
        if method == 'fft':
            wf = propagate_2D_fresnel(wf, propagation_distance*defocus_factor/propagation_steps)
        elif method == 'convolution':
            wf = propagate_2D_fresnel_convolution(wf, propagation_distance*defocus_factor/propagation_steps)
        elif method == 'srw':
            wf = propagate_2D_fresnel_srw(wf, propagation_distance*defocus_factor/propagation_steps)
        elif method == 'fraunhofer':
            wf = propagate_2D_fraunhofer(wf, propagation_distance*defocus_factor/propagation_steps)
        else:
            raise Exception("Not implemented method: %s"%method)




    horizontal_profile = wf.get_intensity()[:,wf.size()[1]/2]
    horizontal_profile /= horizontal_profile.max()
    print("FWHM of the horizontal profile: %g um"%(1e6*line_fwhm(horizontal_profile)*wf.delta()[0]))
    vertical_profile = wf.get_intensity()[wf.size()[0]/2,:]
    vertical_profile /= vertical_profile.max()
    print("FWHM of the vertical profile: %g um"%(1e6*line_fwhm(vertical_profile)*wf.delta()[1]))

    if do_plot:
        from srxraylib.plot.gol import plot,plot_image
        plot_image(wf.get_intensity(),1e6*wf.get_coordinate_x(),1e6*wf.get_coordinate_y(),
                   xtitle="X um",ytitle="Y um",title='intensity (%s)'%method,show=0)
        # plot_image(wf.get_amplitude(),wf.get_coordinate_x(),wf.get_coordinate_y(),title='amplitude (%s)'%method,show=0)
        plot_image(wf.get_phase(),1e6*wf.get_coordinate_x(),1e6*wf.get_coordinate_y(),
                   xtitle="X um",ytitle="Y um",title='phase (%s)'%method,show=0)

        plot(wf.get_coordinate_x(),horizontal_profile,
             wf.get_coordinate_y(),vertical_profile,
             legend=['Horizontal profile','Vertical profile'],title="%s %s"%(plot_title,method),show=show)

    print("Output intensity: ",wf.get_intensity().sum())
    return wf,wf.get_coordinate_x(),horizontal_profile


#
# main function
#

def main(beamline,pixels=100):

    npixels_x = pixels
    npixels_y = npixels_x


    pixelsize_x = beamline['gapH'] / npixels_x
    pixelsize_y = beamline['gapV'] / npixels_y

    print("pixelsize X=%f,Y=%f: "%(pixelsize_x,pixelsize_y))

    propagation_distance = beamline['distance']

    #
    # initialize wavefronts of dimension equal to the lens
    #
    wf_fft = Wavefront2D.initialize_wavefront_from_range(x_min=-pixelsize_x*npixels_x/2,x_max=pixelsize_x*npixels_x/2,
                                                     y_min=-pixelsize_y*npixels_y/2,y_max=pixelsize_y*npixels_y/2,
                                                     number_of_points=(npixels_x,npixels_y),wavelength=1e-10)


    gamma = beamline['ElectronEnergy'] / (codata_mee * 1e-3)
    print ("Gamma: %f \n"%(gamma))


    # photon_wavelength = (1 + beamline['Kv']**2 / 2.0) / 2 / gamma**2 * beamline["PeriodID"] / beamline['harmonicID']

    photon_wavelength = m2ev / beamline["photonEnergy"]


    print ("Photon wavelength [A]: %g \n"%(1e10*photon_wavelength))
    print ("Photon energy [eV]: %g \n"%(beamline["photonEnergy"]))


    myBeam = ElectronBeam(Electron_energy=beamline['ElectronEnergy'], I_current=beamline['ElectronCurrent'])
    myUndulator = MagneticStructureUndulatorPlane(K=beamline['Kv'], period_length=beamline['PeriodID'],
                        length=beamline['PeriodID']*beamline['NPeriods'])


    XX = wf_fft.get_mesh_x()
    YY = wf_fft.get_mesh_y()
    X = wf_fft.get_coordinate_x()
    Y = wf_fft.get_coordinate_y()

    source = SourceUndulatorPlane(undulator=myUndulator,
                        electron_beam=myBeam, magnetic_field=None)


    omega = beamline["photonEnergy"] * codata.e / codata.hbar
    Nb_pts_trajectory = int(source.choose_nb_pts_trajectory(0.01,photon_frequency=omega))
    print("Number of trajectory points: ",Nb_pts_trajectory)


    traj_fact = TrajectoryFactory(Nb_pts=Nb_pts_trajectory,method=TRAJECTORY_METHOD_ODE,
                                  initial_condition=None)

    print("Number of trajectory points: ",traj_fact.Nb_pts)

    if (traj_fact.initial_condition == None):
        traj_fact.initial_condition = source.choose_initial_contidion_automatic()

    print("Number of trajectory points: ",traj_fact.Nb_pts,traj_fact.initial_condition)
    #print('step 2')

    rad_fact = RadiationFactory(method=RADIATION_METHOD_NEAR_FIELD, photon_frequency=omega)


    #print('step 3')
    trajectory = traj_fact.create_from_source(source=source)


    #print('step 4')
    radiation = rad_fact.create_for_one_relativistic_electron(trajectory=trajectory, source=source,
                        XY_are_list=False,distance=beamline['distance'], X=X, Y=Y)

    efield = rad_fact.calculate_electrical_field(trajectory=trajectory,source=source,
                        distance=beamline['distance'],X_array=XX,Y_array=YY)

    tmp = efield.electrical_field()[:,:,0]


    wf_fft.set_photon_energy(beamline["photonEnergy"])

    wf_fft.set_complex_amplitude( tmp )

    # plot

    plot_image(wf_fft.get_intensity(),1e6*wf_fft.get_coordinate_x(),1e6*wf_fft.get_coordinate_y(),
               xtitle="X um",ytitle="Y um",title="UND source at lens plane",show=1)

    # apply lens

    focal_length = propagation_distance / 2

    wf_fft.apply_ideal_lens(focal_length,focal_length)



    plot_image(wf_fft.get_phase(),1e6*wf_fft.get_coordinate_x(),1e6*wf_fft.get_coordinate_y(),
               title="Phase just after the lens",xtitle="X um",ytitle="Y um",show=1)

    wf_fft, x_fft, y_fft = propagation_to_image(wf_fft,do_plot=0,method='fft',
                            propagation_steps=1,
                            propagation_distance = propagation_distance, defocus_factor=1.0)


    plot_image(wf_fft.get_intensity(),1e6*wf_fft.get_coordinate_x(),1e6*wf_fft.get_coordinate_y(),
               title="Intensity at image plane",xtitle="X um",ytitle="Y um",show=1)

    if do_plot:
        plot_table(1e6*x_fft,y_fft,ytitle="Intensity",xtitle="x coordinate [um]",
                   title="Comparison 1:1 focusing ")

if __name__ == "__main__":

    beamline = {}
    # beamline['name'] = "ESRF_NEW_OB"
    # beamline['ElectronBeamDivergenceH'] = 5.2e-6    # these values are not used (zero emittance)
    # beamline['ElectronBeamDivergenceV'] = 1.4e-6    # these values are not used (zero emittance)
    # beamline['ElectronBeamSizeH'] = 27.2e-6         # these values are not used (zero emittance)
    # beamline['ElectronBeamSizeV'] = 3.4e-6          # these values are not used (zero emittance)
    # beamline['ElectronEnergySpread'] = 0.001        # these values are not used (zero emittance)
    beamline['ElectronCurrent'] = 0.5
    beamline['ElectronEnergy'] = 3.0

    beamline['NPeriods'] = 100
    beamline['PeriodID'] = 0.020
    beamline['distance'] =   30.0
    beamline['gapH']      = 0.0016
    beamline['gapV']      = 0.0016


    # get K
    gamma = beamline['ElectronEnergy'] / (codata_mee * 1e-3)
    resonance_energy = 10009.0 / 3
    resonance_wavelength  = m2ev / resonance_energy
    # lambda = ( 1 + K^2/2) lambdau / (2 gamma^2)
    Kv = numpy.sqrt( 2 * ( resonance_wavelength * 2 * gamma**2 / beamline['PeriodID'] - 1) )
    print(">>>>>m2ev=%g Kv=%g  Eres=%g  lambda_reson=%g"%(m2ev,Kv,resonance_energy,resonance_wavelength))

    beamline['Kv'] = Kv
    beamline["photonEnergy"] = 10030.0 # 9996.0 # 10009.0


    main(beamline,pixels=300)