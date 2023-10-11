
def get_source(do_plot=0):
    #
    # script to make the calculations (created by XOPPY:wiggler_radiation)
    #

    from xoppylib.sources.xoppy_bm_wiggler import xoppy_calc_wiggler_radiation

    h5_parameters = dict()
    h5_parameters["ELECTRONENERGY"]          = 6.0
    h5_parameters["ELECTRONCURRENT"]         = 0.2
    h5_parameters["PERIODID"]                = 0.12
    h5_parameters["NPERIODS"]                = 37.0
    h5_parameters["KV"]                      = 22.416
    h5_parameters["FIELD"]                   = 1   # 0= sinusoidal, 1=from file
    h5_parameters["FILE"]                    = '/nobackup/gurb1/srio/Oasys/ShadowOui-Tutorial/SOS-WORKSHOP/EBS-WIGGLERS/SW_BM18.txt'
    h5_parameters["POLARIZATION"]            = 0  # 0=total, 1=s, 2=p
    h5_parameters["DISTANCE"]                = 30.0
    h5_parameters["HSLITPOINTS"]             = 50
    h5_parameters["VSLITPOINTS"]             = 50
    h5_parameters["PHOTONENERGYMIN"]         = 100.0
    h5_parameters["PHOTONENERGYMAX"]         = 200000.0
    h5_parameters["PHOTONENERGYPOINTS"]      = 500
    h5_parameters["SHIFT_X_FLAG"]            = 5
    h5_parameters["SHIFT_X_VALUE"]           = 0.0
    h5_parameters["SHIFT_BETAX_FLAG"]        = 4
    h5_parameters["SHIFT_BETAX_VALUE"]       = 0.0
    h5_parameters["CONVOLUTION"]             = 1
    h5_parameters["PASSEPARTOUT"]            = 3.0

    energy, horizontal, vertical, flux3D, traj = xoppy_calc_wiggler_radiation(
            ELECTRONENERGY           = h5_parameters["ELECTRONENERGY"]         ,
            ELECTRONCURRENT          = h5_parameters["ELECTRONCURRENT"]        ,
            PERIODID                 = h5_parameters["PERIODID"]               ,
            NPERIODS                 = h5_parameters["NPERIODS"]               ,
            KV                       = h5_parameters["KV"]                     ,
            FIELD                    = h5_parameters["FIELD"]                  ,
            FILE                     = h5_parameters["FILE"]                   ,
            POLARIZATION             = h5_parameters["POLARIZATION"]           ,
            DISTANCE                 = h5_parameters["DISTANCE"]               ,
            HSLITPOINTS              = h5_parameters["HSLITPOINTS"]            ,
            VSLITPOINTS              = h5_parameters["VSLITPOINTS"]            ,
            PHOTONENERGYMIN          = h5_parameters["PHOTONENERGYMIN"]        ,
            PHOTONENERGYMAX          = h5_parameters["PHOTONENERGYMAX"]        ,
            PHOTONENERGYPOINTS       = h5_parameters["PHOTONENERGYPOINTS"]     ,
            SHIFT_X_FLAG             = h5_parameters["SHIFT_X_FLAG"]           ,
            SHIFT_X_VALUE            = h5_parameters["SHIFT_X_VALUE"]          ,
            SHIFT_BETAX_FLAG         = h5_parameters["SHIFT_BETAX_FLAG"]       ,
            SHIFT_BETAX_VALUE        = h5_parameters["SHIFT_BETAX_VALUE"]      ,
            CONVOLUTION              = h5_parameters["CONVOLUTION"]            ,
            PASSEPARTOUT             = h5_parameters["PASSEPARTOUT"]            ,
            h5_file                  = "wiggler_radiation.h5"                  ,
            h5_entry_name            = "XOPPY_RADIATION"                       ,
            h5_initialize            = True                                    ,
            h5_parameters            = h5_parameters                           ,
            )

    if do_plot:
        # example plot
        from srxraylib.plot.gol import plot_image
        plot_image(flux3D[0],horizontal,vertical,title="Flux [photons/s] per 0.1 bw per mm2 at %9.3f eV"%(100.0),xtitle="H [mm]",ytitle="V [mm]")
    #
    # end script
    #
    return energy, horizontal, vertical, flux3D

def get_source_from_file(do_plot=0):
    #
    # script to load the wiggler radiation calculations from file wiggler_radiation.h5
    #
    import h5py
    hf = h5py.File('/tmp_14_days/srio/wiggler_radiation.h5', 'r')

    flux3D = hf["/XOPPY_RADIATION/Radiation/stack_data"][:]
    energy = hf["/XOPPY_RADIATION/Radiation/axis0"][:]
    horizontal = hf["/XOPPY_RADIATION/Radiation/axis1"][:]
    vertical = hf["/XOPPY_RADIATION/Radiation/axis2"][:]
    traj = hf["/XOPPY_RADIATION/trajectory/traj"][:]

    hf.close()

    # example plot
    if do_plot:
        from srxraylib.plot.gol import plot_image
        plot_image(flux3D[0], horizontal, vertical, title="Flux [photons/s] per 0.1 bw per mm2 at %9.3f eV" % (100.0),
                   xtitle="H [mm]", ytitle="V [mm]")
    #
    # end script
    #
    return energy, horizontal, vertical, flux3D

def slit(energy, horizontal, vertical, flux3D, do_plot=0):
    #
    # script to make the calculations (created by XOPPY:power3Dcomponent)
    #

    import numpy
    from xoppylib.power.power3d import calculate_component_absorbance_and_transmittance
    from xoppylib.power.power3d import apply_transmittance_to_incident_beam

    # compute local transmittance and absorbance
    e0, h0, v0, f0  = energy, horizontal, vertical, flux3D
    transmittance, absorbance, E, H, V, txt = calculate_component_absorbance_and_transmittance(
                    e0, # energy in eV
                    h0, # h in mm
                    v0, # v in mm
                    substance='Be',
                    thick=0.5,
                    angle=3.0,
                    defection=1,
                    dens='?',
                    roughness=0.0,
                    flags=2,
                    hgap=70.0,
                    vgap=6.0,
                    hgapcenter=0.0,
                    vgapcenter=0.0,
                    hmag=1.0,
                    vmag=1.0,
                    hrot=0.0,
                    vrot=0.0,
                    )

    # apply transmittance to incident beam
    f_transmitted, e, h, v = apply_transmittance_to_incident_beam(transmittance, f0, e0, h0, v0,
                                      flags = 2,
                                      hgap = 70.0,
                                      vgap = 6.0,
                                      hgapcenter = 0.0,
                                      vgapcenter = 0.0,
                                      hmag = 1.0,
                                      vmag = 1.0,
                                      interpolation_flag     = 0,
                                      interpolation_factor_h = 1.0,
                                      interpolation_factor_v = 1.0,
                                      slit_crop = 0,
                                    )

    f_absorbed = f0 * absorbance / (H[0] / h0[0]) / (V[0] / v0[0])

    # data to pass
    energy, horizontal, vertical, flux3D = e, h, v, f_transmitted

    #
    # example plots
    #
    if do_plot:
        from srxraylib.plot.gol import plot_image
        import scipy.constants as codata
        from xoppylib.power.power3d import integral_2d

        # transmitted/reflected beam

        spectral_power_transmitted = f_transmitted * codata.e * 1e3
        plot_image(spectral_power_transmitted[0,:,:],h,v,title="Transmitted Spectral Power Density [W/eV/mm2] at E=%g eV" % (100.0),xtitle="H [mm]",ytitle="V [mm]",aspect='auto')

        power_density_transmitted = numpy.trapz(spectral_power_transmitted, e, axis=0)
        power_density_integral = integral_2d(power_density_transmitted, h, v)
        plot_image(power_density_transmitted, h, v,
                         xtitle='H [mm] (normal to beam)',
                         ytitle='V [mm] (normal to beam)',
                         title='Power Density [W/mm^2]. Integral: %6.3f W'%power_density_integral,aspect='auto')

        # local absorption

        spectral_power_density_absorbed = f_absorbed * codata.e * 1e3

        plot_image(spectral_power_density_absorbed[0,:,:],H,V,title="Absorbed Spectral Power Density [W/eV/mm2] at E=%g eV" % (100.0),xtitle="H [mm]",ytitle="V [mm]",aspect='auto')

        power_density_absorbed = numpy.trapz(spectral_power_density_absorbed, E, axis=0)
        power_density_integral = integral_2d(power_density_absorbed, H, V)
        plot_image(power_density_absorbed, H, V,
                         xtitle='H [mm] (o.e. coordinates)',
                         ytitle='V [mm] (o.e. coordinates)',
                         title='Absorbed Power Density [W/mm^2]. Integral: %6.3f W'%power_density_integral,aspect='auto')

    #
    # end script
    #
    return energy, horizontal, vertical, flux3D

def filter(energy, horizontal, vertical, flux3D, do_plot=0,
           substance='C',
           thick=0.2,
           dens='?',
           ):
    #
    # script to make the calculations (created by XOPPY:power3Dcomponent)
    #

    import numpy
    from xoppylib.power.power3d import calculate_component_absorbance_and_transmittance
    from xoppylib.power.power3d import apply_transmittance_to_incident_beam

    # compute local transmittance and absorbance
    e0, h0, v0, f0  = energy, horizontal, vertical, flux3D
    transmittance, absorbance, E, H, V, txt = calculate_component_absorbance_and_transmittance(
                    e0, # energy in eV
                    h0, # h in mm
                    v0, # v in mm
                    substance=substance, #****
                    thick=thick,         #****
                    angle=3.0,
                    defection=1,
                    dens=dens,
                    roughness=0.0,
                    flags=0,
                    hgap=1000.0,
                    vgap=1000.0,
                    hgapcenter=0.0,
                    vgapcenter=0.0,
                    hmag=1.0,
                    vmag=1.0,
                    hrot=0.0,
                    vrot=0.0,
                    )

    # apply transmittance to incident beam
    f_transmitted, e, h, v = apply_transmittance_to_incident_beam(transmittance, f0, e0, h0, v0,
                                      flags = 0,
                                      hgap = 1000.0,
                                      vgap = 1000.0,
                                      hgapcenter = 0.0,
                                      vgapcenter = 0.0,
                                      hmag = 1.0,
                                      vmag = 1.0,
                                      interpolation_flag     = 0,
                                      interpolation_factor_h = 1.0,
                                      interpolation_factor_v = 1.0,
                                      slit_crop = 0,
                                    )

    f_absorbed = f0 * absorbance / (H[0] / h0[0]) / (V[0] / v0[0])

    # data to pass
    energy, horizontal, vertical, flux3D = e, h, v, f_transmitted

    #
    # example plots
    #
    if do_plot:
        from srxraylib.plot.gol import plot_image
        import scipy.constants as codata
        from xoppylib.power.power3d import integral_2d

        # transmitted/reflected beam

        spectral_power_transmitted = f_transmitted * codata.e * 1e3
        plot_image(spectral_power_transmitted[0,:,:],h,v,title="Transmitted Spectral Power Density [W/eV/mm2] at E=%g eV" % (100.0),xtitle="H [mm]",ytitle="V [mm]",aspect='auto')

        power_density_transmitted = numpy.trapz(spectral_power_transmitted, e, axis=0)
        power_density_integral = integral_2d(power_density_transmitted, h, v)
        plot_image(power_density_transmitted, h, v,
                         xtitle='H [mm] (normal to beam)',
                         ytitle='V [mm] (normal to beam)',
                         title='Power Density [W/mm^2]. Integral: %6.3f W'%power_density_integral,aspect='auto')

        # local absorption

        spectral_power_density_absorbed = f_absorbed * codata.e * 1e3

        plot_image(spectral_power_density_absorbed[0,:,:],H,V,title="Absorbed Spectral Power Density [W/eV/mm2] at E=%g eV" % (100.0),xtitle="H [mm]",ytitle="V [mm]",aspect='auto')

        power_density_absorbed = numpy.trapz(spectral_power_density_absorbed, E, axis=0)
        power_density_integral = integral_2d(power_density_absorbed, H, V)
        plot_image(power_density_absorbed, H, V,
                         xtitle='H [mm] (o.e. coordinates)',
                         ytitle='V [mm] (o.e. coordinates)',
                         title='Absorbed Power Density [W/mm^2]. Integral: %6.3f W'%power_density_integral,aspect='auto')

    #
    # end script
    #
    return energy, horizontal, vertical, flux3D

def magnifier(energy, horizontal, vertical, flux3D, do_plot=0):
    #
    # script to make the calculations (created by XOPPY:power3Dcomponent)
    #

    import numpy
    from xoppylib.power.power3d import calculate_component_absorbance_and_transmittance
    from xoppylib.power.power3d import apply_transmittance_to_incident_beam

    # compute local transmittance and absorbance
    e0, h0, v0, f0  = energy, horizontal, vertical, flux3D
    transmittance, absorbance, E, H, V, txt = calculate_component_absorbance_and_transmittance(
                    e0, # energy in eV
                    h0, # h in mm
                    v0, # v in mm
                    substance='C',
                    thick=0.2,
                    angle=3.0,
                    defection=1,
                    dens='?',
                    roughness=0.0,
                    flags=3,
                    hgap=1000.0,
                    vgap=1000.0,
                    hgapcenter=0.0,
                    vgapcenter=0.0,
                    hmag=172.055 / 26.276,
                    vmag=172.055 / 26.276,
                    hrot=0.0,
                    vrot=0.0,
                    )

    # apply transmittance to incident beam
    f_transmitted, e, h, v = apply_transmittance_to_incident_beam(transmittance, f0, e0, h0, v0,
                                      flags = 3,
                                      hgap = 1000.0,
                                      vgap = 1000.0,
                                      hgapcenter = 0.0,
                                      vgapcenter = 0.0,
                                      hmag = 6.54,
                                      vmag = 6.54,
                                      interpolation_flag     = 0,
                                      interpolation_factor_h = 1.0,
                                      interpolation_factor_v = 1.0,
                                      slit_crop = 0,
                                    )

    f_absorbed = f0 * absorbance / (H[0] / h0[0]) / (V[0] / v0[0])

    # data to pass
    energy, horizontal, vertical, flux3D = e, h, v, f_transmitted

    #
    # example plots
    #
    if do_plot:
        from srxraylib.plot.gol import plot_image
        import scipy.constants as codata
        from xoppylib.power.power3d import integral_2d

        # transmitted/reflected beam

        spectral_power_transmitted = f_transmitted * codata.e * 1e3
        plot_image(spectral_power_transmitted[0,:,:],h,v,title="Transmitted Spectral Power Density [W/eV/mm2] at E=%g eV" % (100.0),xtitle="H [mm]",ytitle="V [mm]",aspect='auto')

        power_density_transmitted = numpy.trapz(spectral_power_transmitted, e, axis=0)
        power_density_integral = integral_2d(power_density_transmitted, h, v)
        plot_image(power_density_transmitted, h, v,
                         xtitle='H [mm] (normal to beam)',
                         ytitle='V [mm] (normal to beam)',
                         title='Power Density [W/mm^2]. Integral: %6.3f W'%power_density_integral,aspect='auto')

        # local absorption

        spectral_power_density_absorbed = f_absorbed * codata.e * 1e3

        plot_image(spectral_power_density_absorbed[0,:,:],H,V,title="Absorbed Spectral Power Density [W/eV/mm2] at E=%g eV" % (100.0),xtitle="H [mm]",ytitle="V [mm]",aspect='auto')

        power_density_absorbed = numpy.trapz(spectral_power_density_absorbed, E, axis=0)
        power_density_integral = integral_2d(power_density_absorbed, H, V)
        plot_image(power_density_absorbed, H, V,
                         xtitle='H [mm] (o.e. coordinates)',
                         ytitle='V [mm] (o.e. coordinates)',
                         title='Absorbed Power Density [W/mm^2]. Integral: %6.3f W'%power_density_integral,aspect='auto')

    #
    # end script
    #
    return energy, horizontal, vertical, flux3D

if __name__ == "__main__":

    e, h, v, f = get_source_from_file()
    e, h, v, f = slit  (e, h, v, f)

    filter_substance = ['Ti23.95Al2.83V','C2H6OSi','SiO2','C2H6OSi','Ti23.95Al2.83V', ]
    filter_thick     = [6.51,            13.0,     10.0,  13.0,     6.51, ]
    filter_dens      = [4.43,            0.9493,   2.32,  0.9493,   4.43, ]
    for i in range(len(filter_substance)):
        e, h, v, f = filter(e, h, v, f,
                            substance=filter_substance[i],
                            thick=filter_thick[i],
                            dens=filter_dens[i])

    e, h, v, f_transmitted = magnifier(e, h, v, f, do_plot=0)


    import numpy
    import scipy.constants as codata
    from srxraylib.plot.gol import plot_image, plot
    from xoppylib.power.power3d import integral_2d

    spectral_power_transmitted = f_transmitted * codata.e * 1e3
    power_density_transmitted = numpy.trapz(spectral_power_transmitted, e, axis=0)
    power_density_integral = integral_2d(power_density_transmitted, h, v)
    plot_image(power_density_transmitted, h, v,
               xtitle='H [mm] (normal to beam)',
               ytitle='V [mm] (normal to beam)',
               title='Power Density [W/mm^2]. Integral: %6.3f W' % power_density_integral, aspect='auto', show=0)


    f_transmitted_per_cte_bw = numpy.zeros_like(f_transmitted)
    for i in range(f_transmitted.shape[0]):
        f_transmitted_per_cte_bw[i, :, :] = f_transmitted[i, :, :] / (1e-3 * e[i])

    flux_density_transmitted = numpy.trapz(f_transmitted_per_cte_bw, e, axis=0)

    flux_density_integral = integral_2d(flux_density_transmitted, h, v)
    plot_image(flux_density_transmitted, h, v,
               xtitle='H [mm] (normal to beam)',
               ytitle='V [mm] (normal to beam)',
               title='Flux Density [ph/s/mm^2]. Integral: %g photons/s' % flux_density_integral, aspect='auto', show=0)


    spectrum = numpy.trapz(numpy.trapz(f_transmitted, v, axis=2), h, axis=1)
    plot(e, spectrum, xlog=0, ylog=0, ytitle="Photons/s/0.1%bw", show=0)

    spectrum_eV = numpy.trapz(numpy.trapz(f_transmitted_per_cte_bw, v, axis=2), h, axis=1)
    plot(e, spectrum_eV, xlog=0, ylog=0, ytitle="Photons/s/eV")