import numpy
import os
import Shadow

import sys
sys.path.append('/scisoft/xop2.4/extensions/shadowvui/python_scripts/')
from transfocator_id30b_lighted import transfocator_compute_configuration
from srundplug import getBeamline,calc1dUrgent,calc1dSrw

def run_source(label='ESRF_LB_OB',undulator='U20',photon_energy_ev=14200.0,\
        photon_energy_bandwidth=1.0,compute_flux=1):

    #myBL = getBeamline('ESRF_LB_OB')
    #myBL = getBeamline('ESRF_LB')
    #myBL = getBeamline('ESRF_NEW_OB')

    myBL = getBeamline(label)

    SIGMAX=myBL['ElectronBeamSizeH']
    SIGMAZ= myBL['ElectronBeamSizeV']
    SIGDIX=  myBL['ElectronBeamDivergenceH']
    SIGDIZ= myBL['ElectronBeamDivergenceV'] 

    if undulator == 'U20':
        myBL['Kv'] = 0.63
        myBL['PeriodID'] = 0.0202 
        myBL['N'] = 70 

    if undulator == 'U14':
        myBL['Kv'] = 1.2
        myBL['PeriodID'] = 0.014
        myBL['N'] = 143

    undulator_length_in_m= myBL['PeriodID'] * myBL['N'] 
    print("UNDULATOR LENGTH: ",undulator_length_in_m)

    #SIGMAX=1e-5
    #SIGMAZ= 1e-5
    #SIGDIX=  1e-6
    #SIGDIZ= 1e-6

    #
    #SHADOW SOURCE
    #

    # set Gaussian source
    src = Shadow.Source()
    src.set_energy_monochromatic(photon_energy_ev)
    src.set_gauss(SIGMAX*1e2,SIGMAZ*1e2,SIGDIX,SIGDIZ)

    print("\n\nElectron sizes stored H:%f um, V:%f um;\nelectron divergences: H:%f urad, V:%f urad"%\
          (src.SIGMAX*1e4,src.SIGMAZ*1e4,src.SIGDIX*1e6,src.SIGDIZ*1e6))

    src.apply_gaussian_undulator(undulator_length_in_m=undulator_length_in_m, user_unit_to_m=1e-2, verbose=1)

    print("\n\nElectron sizes stored (undulator) H:%f um, V:%f um;\nelectron divergences: H:%f urad, V:%f urad"%\
          (src.SIGMAX*1e4,src.SIGMAZ*1e4,src.SIGDIX*1e6,src.SIGDIZ*1e6))

    src.set_energy_box(photon_energy_ev-0.5*photon_energy_bandwidth,photon_energy_ev+0.5*photon_energy_bandwidth)

    print("SIGMAX = %g"%src.SIGMAX)
    print("SIGMAZ = %g"%src.SIGMAZ)
    print("SIGDIX = %g"%src.SIGDIX)
    print("SIGDIZ = %g"%src.SIGDIZ)

    if compute_flux:
        myBL['d'] = 28.20
        myBL['gapH'] =  myBL['d'] * 6.0 * src.SIGDIX
        myBL['gapV'] =  myBL['d'] * 6.0 * src.SIGDIZ
        myBL['ElectronEnergySpread'] = 0.001
        eUrgent,fUrgent = calc1dUrgent(myBL, \
            photonEnergyMin=10000.0,photonEnergyMax=18000.0,
            photonEnergyPoints=250,fileName=undulator+'.spec',fileAppend=False)
        #eSrw,fSrw = calc1dSrw(myBL, \
        #    photonEnergyMin=10000.0,photonEnergyMax=18000.0,
        #    photonEnergyPoints=250,fileName=undulator+'.spec',fileAppend=True)
        flux_1ev_bandwidth = numpy.max(fUrgent) / (photon_energy_ev*1e-3)
    else:
        flux_1ev_bandwidth = 1.0

    src.NPOINT = 250000
    src.ISTAR1 = 110011 

    src.write("start.00")


    # create source
    beam = Shadow.Beam()
    beam.genSource(src)
    beam.write("begin.dat")
    src.write("end.00")

    return (beam,src,flux_1ev_bandwidth)

#
# KB beamline
#
def set_beamline_kb(configuration=0):
    # crystal (from source to exit slit plane)
    if mono:
        oe1 = Shadow.OE()
        oe1.F_CRYSTAL = 1
        oe1.ALPHA = 90.0 
        oe1.T_SOURCE = 3000.0
        oe1.T_IMAGE = 1090.0
        oe1.FMIRR = 5
        oe1.F_CENTRAL = 1
        oe1.F_PHOT_CENT = 0 #eV
        oe1.PHOT_CENT = photon_energy_ev
        oe1.FILE_REFL = "Si5_55.111".encode('utf-8')
        oe1.FWRITE = 0
    
        # screen before mono (screen=1, index=0)
        if 0:
            oe1.F_SCREEN = 1
            oe1.N_SCREEN = 1
            oe1.I_SCREEN[0] = 1 # 1=Before, 0=After
            oe1.SL_DIS[0] = 180.0
            oe1.FWRITE = 0 # write all files 

    # KB
    if mono:
        dist_mono = 0
        m_o_a = 90.0
    else:
        dist_mono = 4090.0
        m_o_a = 0.0

    kb_separation = 14.0
    d_sample = 4515.0
    d_slit = 4090.0

    if configuration == 6:
        # KB center at 4315
        d_from_sample = 200.0
    if configuration == 7:
        # KB center at 4515-85 = 4430 from sample
        d_from_sample = 85.0

    d_center = d_sample - d_from_sample

    kb = Shadow.CompoundOE(name='KB')
 
    kb.append_kb(dist_mono+(d_center-d_slit),(d_sample-d_center),\
        separation=kb_separation,\
        mirror_orientation_angle=m_o_a,\
        grazing_angles_mrad=[3.9,17.8],\
        focal_positions = [d_center,d_sample-d_center],\
        shape=[2,2],dimensions1=[6,24],dimensions2=[6,14],\
        reflectivity_kind=[0,0],reflectivity_files=["",""],\
        surface_error_files=["waviness_0p4urad.dat","waviness_0p4urad.dat"])

    print("    d_kb_from_sample: %6.2f cm, averaged demagnification: %.2f "%( d_from_sample,(d_sample-d_from_sample)/d_from_sample)) 
    print("    inter-mirror distance: %6.2f cm "%( kb_separation)) 
    if ihit: itmp = input("Hit a key to continue...")

    #
    # build beamline
    #
    bl = Shadow.CompoundOE(name='ID23-2')
    if mono:
        bl.append(oe1)
    bl.append(kb)
    bl.append(Shadow.OE().set_empty(ALPHA=-90))

    return bl


#
# 2D crl beamline
#

def set_beamline_twoCRLs(configuration=0):
    # crystal (from source to exit slit plane)
    if mono:
        oe1 = Shadow.OE()
        oe1.F_CRYSTAL = 1
        oe1.ALPHA = 90.0 
        oe1.T_SOURCE = 3000.0
        oe1.T_IMAGE = 1090.0
        oe1.FMIRR = 5
        oe1.F_CENTRAL = 1
        oe1.F_PHOT_CENT = 0 #eV
        oe1.PHOT_CENT = photon_energy_ev
        oe1.FILE_REFL = "Si5_55.111".encode('utf-8')
        oe1.FWRITE = 0
    
        # screen before mono (screen=1, index=0)
        if 0:
            oe1.F_SCREEN = 1
            oe1.N_SCREEN = 1
            oe1.I_SCREEN[0] = 1 # 1=Before, 0=After
            oe1.SL_DIS[0] = 180.0
            oe1.FWRITE = 0 # write all files 

    d_sample = 4515.0 
    d_slit = 4090.0

    # TF
    if configuration == 1:   # Be R=200um 
        d_crlv_from_sample = 1475.0
        d_crlh_from_sample = 251.0
    if configuration == 2:         # Be R=200um 
        d_crlv_from_sample = 302.0 
        d_crlh_from_sample = 100.9 
    if configuration == 3:         # Be R=200um 
        d_crlv_from_sample = 302.0 
        d_crlh_from_sample = 100.9 

    d_continuation_plane = d_sample - 0.5 * (d_crlv_from_sample + d_crlh_from_sample)

    #vertical foci
    d_centerV = d_sample - d_crlv_from_sample 
    f_target = 1.0/(1.0/d_crlv_from_sample + 1.0/d_centerV)
    delta_be_14200ev = 1.68692763813e-06
    nlensesV = (200e-4) / 2 /delta_be_14200ev / f_target
    print("V -> q:%.2f, q:%.2f, f:%.2f, Nv:%.3f \n"%(d_crlv_from_sample,d_centerV,f_target,nlensesV))
    #horizontal foci
    d_centerH = d_sample - d_crlh_from_sample 
    f_target = 1.0/(1.0/d_crlh_from_sample + 1.0/d_centerH)
    delta_be_14200ev = 1.68692763813e-06
    nlensesH = (200e-4) / 2 /delta_be_14200ev / f_target
    print("H -> q:%.2f, p:%.2f, f:%.2f, Nh:%.3f \n"%(d_crlh_from_sample,d_centerH,f_target,nlensesH))

    nlensesV = int(numpy.round(nlensesV))
    nlensesH = int(numpy.round(nlensesH))

    print("    d_crlv_from_sample: %6.2f cm, V demagnification: %.2f, N lenses V:%d "%( d_crlv_from_sample,d_centerV/d_crlv_from_sample,nlensesV)) 
    print("    d_crlh_from_sample: %6.2f cm, H demagnification: %.2f, N lenses H:%d "%( d_crlh_from_sample,d_centerH/d_crlh_from_sample,nlensesV)) 

    if ihit: itmp = input("Hit a key to continue...")


    # vertical CRL
    crlV = Shadow.CompoundOE(name='CRLv')
    crlV.append_crl(0.0,0.0, nlenses=nlensesV, radius=200e-4, \
            thickness=1000e-4, interthickness=30e-4, \
            diameter=[2500e-4,894e-4], \
            surface_shape=4, convex_to_the_beam=0, cylinder_angle=0,\
            prerefl_file="Be2_55.dat", use_ccc=0)
    crlV_length = crlV.length()

    # horizontal CRL
    crlH = Shadow.CompoundOE(name='CRLh')
    crlH.append_crl(0.0,0.0, nlenses=nlensesH, radius=200e-4, \
            thickness=1000e-4, interthickness=30e-4, \
            diameter=[2500e-4,894e-4], \
            surface_shape=4, convex_to_the_beam=0, cylinder_angle=90,\
            prerefl_file="Be2_55.dat", use_ccc=0)
    crlH_length = crlH.length()


    print("CRLv length: %f cm"%(crlV_length))
    print("CRLh length: %f cm"%(crlH_length))
    #itmp = input("Hit a key to continue...")

    if mono:
        crlV.add_drift_space_upstream(d_centerV-d_slit-0.5*crlV_length)
    else:
        crlV.add_drift_space_upstream(d_slit+d_centerV-d_slit-0.5*crlV_length)
    crlV.add_drift_space_downstream(d_continuation_plane-d_centerV-0.5*crlV_length)

    crlH.add_drift_space_upstream(d_centerH-d_continuation_plane-0.5*crlH_length)
    crlH.add_drift_space_downstream(d_sample-d_centerH-0.5*crlH_length)

    #
    # build beamline
    #

    bl = Shadow.CompoundOE(name='ID23-2')
    if mono:
        bl.append(oe1)
        oe0 = Shadow.OE()
        oe0.set_empty(ALPHA=-90)
        bl.append(oe0) # switch axes

    bl.append(crlV)
    bl.append(crlH)
    print(bl.info())
    #itmp = input("Hit a key to continue...")

    return bl

#
# 2D crl beamline
#

def set_beamline_crl2D(configuration=0):
    # crystal (from source to exit slit plane)
    if mono:
        oe1 = Shadow.OE()
        oe1.F_CRYSTAL = 1
        oe1.ALPHA = 90.0 
        oe1.T_SOURCE = 3000.0
        oe1.T_IMAGE = 1090.0
        oe1.FMIRR = 5
        oe1.F_CENTRAL = 1
        oe1.F_PHOT_CENT = 0 #eV
        oe1.PHOT_CENT = photon_energy_ev
        oe1.FILE_REFL = "Si5_55.111".encode('utf-8')
        oe1.FWRITE = 0
    
        # screen before mono (screen=1, index=0)
        if 0:
            oe1.F_SCREEN = 1
            oe1.N_SCREEN = 1
            oe1.I_SCREEN[0] = 1 # 1=Before, 0=After
            oe1.SL_DIS[0] = 180.0
            oe1.FWRITE = 0 # write all files 

    d_sample = 4515.0 
    d_slit = 4090.0

    # TF
    if configuration == 1:   # Be R=200um 
        d_crlv_from_sample = 49.95 
    if configuration == 2:         # Be R=200um 
        d_crlv_from_sample = 200.0 
    if configuration == 3:         # Be R=200um 
        d_crlv_from_sample = 302.0


    d_center = d_sample - d_crlv_from_sample 
    f_target = 1.0/(1.0/d_crlv_from_sample + 1.0/d_center)
    delta_be_14200ev = 1.68692763813e-06
    nlenses = (200e-4) / 2 /delta_be_14200ev / f_target
    print("p:%f.2, q:%.2f, f:%.2f, N:%.3f \n"%(d_crlv_from_sample,d_center,f_target,nlenses))
    nlenses = int(numpy.round(nlenses))
    #itmp = input("Hit a key to continue...")


    crl = Shadow.CompoundOE(name='CRL')
    crl.append_crl(0.0,0.0, nlenses=nlenses, radius=200e-4, \
            thickness=1000e-4, interthickness=30e-4, \
            diameter=894e-4, \
            surface_shape=4, convex_to_the_beam=0, cylinder_angle=None,\
            prerefl_file="Be2_55.dat", use_ccc=0)
    crl_length = crl.length()
    print("CRL length: %f cm"%(crl_length))
    if mono:
        crl.add_drift_space_upstream(d_center-d_slit-0.5*crl_length)
    else:
        crl.add_drift_space_upstream(d_slit+d_center-d_slit-0.5*crl_length)
    crl.add_drift_space_downstream(d_crlv_from_sample-0.5*crl_length)

    #
    # build beamline
    #

    bl = Shadow.CompoundOE(name='ID23-2')
    if mono:
        bl.append(oe1)
        oe0 = Shadow.OE()
        oe0.set_empty(ALPHA=-90)
        bl.append(oe0) # switch axes

    bl.append(crl)
    print(bl.info())
    #itmp = input("Hit a key to continue...")

    return bl


def set_beamline_ml(configuration=0):
    # crystal (from source to exit slit plane)
    if mono:
        oe1 = Shadow.OE()
        oe1.F_CRYSTAL = 1
        oe1.ALPHA = 90.0 
        oe1.T_SOURCE = 3000.0
        oe1.T_IMAGE = 1090.0
        oe1.FMIRR = 5
        oe1.F_CENTRAL = 1
        oe1.F_PHOT_CENT = 0 #eV
        oe1.PHOT_CENT = photon_energy_ev
        oe1.FILE_REFL = "Si5_55.111".encode('utf-8')
        oe1.FWRITE = 0
    
        # screen before mono (screen=1, index=0)
        if 0:
            oe1.F_SCREEN = 1
            oe1.N_SCREEN = 1
            oe1.I_SCREEN[0] = 1 # 1=Before, 0=After
            oe1.SL_DIS[0] = 180.0
            oe1.FWRITE = 0 # write all files 

    d_sample = 4515.0 
    d_slit = 4090.0
    if configuration == 4:   # Be R=200um 
        d_crl_from_sample = 1475 # 358 # 251.0
        d_ml_from_sample = 250 # 100.0

    if configuration == 5:   # Be R=200um 
        d_crl_from_sample = 515 # 358 # 251.0
        d_ml_from_sample = 35 # 100.0

    #
    # CRL
    #
    d_center = d_sample - d_crl_from_sample
    f_target = 1.0/(1.0/d_crl_from_sample + 1.0/d_center)
    delta_be_14200ev = 1.68692763813e-06
    nlenses = (200e-4) / 2 /delta_be_14200ev / f_target
    print("p:%f.2, q:%.2f, f:%.2f, N:%.3f \n"%(d_center,d_crl_from_sample,f_target,nlenses))
    print("nlenses: %f"%(nlenses))
    nlenses = int(numpy.round(nlenses))


    crl = Shadow.CompoundOE(name='CRL')
    crl.append_crl(0.0,0.0, nlenses=nlenses, radius=200e-4, \
            thickness=1000e-4, interthickness=30e-4, \
            diameter=[2500e-4,894e-4], \
            surface_shape=4, convex_to_the_beam=0, cylinder_angle=0,\
            prerefl_file="Be2_55.dat", use_ccc=0)
    crl_length = crl.length()
    print("CRL length: %f cm"%(crl_length))
    crl.add_drift_space_upstream(d_center-d_slit-0.5*crl_length)
    crl.add_drift_space_downstream( (d_crl_from_sample-d_ml_from_sample)\
                                     -0.5*crl_length)

    #
    #multilayer
    #
    inc_angle_deg = 90.0 - 15e-3*180/numpy.pi
    d_center = d_sample - d_ml_from_sample

    ml = Shadow.OE()
    ml.FMIRR = 2 # elliptical
    ml.FCYL = 1  # cylindrical
    ml.CIL_ANG = 0.0 # tangential focusing
    #internal parameters
    ml.F_EXT = 0     # internal/calculated
    ml.F_DEFAULT = 0 # focii not coincident with continuation plane
    ml.SSOUR = d_center
    ml.SIMAG = d_ml_from_sample
    ml.THETA = inc_angle_deg
    #dimensions
    ml.FHIT_C = 1    # mirror dimensions finite: yes (1), no(0).
    ml.FSHAPE = 1    # rectangle
    ml.RWIDX1 = 0.5 * 5
    ml.RWIDX2 = 0.5 * 5
    ml.RLEN1  = 0.5 * 24
    ml.RLEN2  = 0.5 * 24
    #distances and orientation angles
    ml.T_SOURCE = 0.0
    ml.T_IMAGE  = 0.0
    ml.T_INCIDENCE = inc_angle_deg
    ml.T_REFLECTION = inc_angle_deg
    ml.ALPHA = 90.0
    # slope errors
    ml.F_RIPPLE = 1 # 0=No
    ml.FILE_RIP = "waviness_0p4urad.dat".encode('utf-8')
    ml.F_G_S = 2
    #set to write all files (overwritten later?)	
    ml.FWRITE = 0

    print("    d_crl_from_sample: %6.2f cm, V demagnification: %.2f, N lenses V:%d "%( d_crl_from_sample,(d_sample-d_crl_from_sample)/d_crl_from_sample,nlenses)) 
    print("    d_ml_from_sample: %6.2f cm, H demagnification: %.2f "%( d_ml_from_sample,(d_sample-d_ml_from_sample)/d_ml_from_sample)) 
    if ihit: itmp = input("Hit a key to continue...")
    #itmp = input("Hit a key to continue...")

    # empty element to rotate axes
    oe0 = Shadow.OE()
    oe0.set_empty(ALPHA=-90)
    #
    # build beamline
    #
    bl = Shadow.CompoundOE(name='ID23-2')
    if mono:
        bl.append(oe1)
        bl.append(oe0) # switch axes
    else:
        crl.add_drift_space_upstream(d_slit)

    bl.append(crl)
    bl.append(ml)

    bl.append(oe0) # switch axes
    bl.add_drift_space_downstream(d_ml_from_sample) # warning: modifies oe0!!

    return bl

def main():
   

    (beam,src,flux_1ev_bandwidth) = run_source(ring,undulator=undulator,\
        photon_energy_ev=photon_energy_ev,\
        photon_energy_bandwidth=photon_energy_bandwidth)

    tk_s_h = beam.histo1(col=1, nbins = 150, nolost=1, ref=1)
    tk_s_v = beam.histo1(col=3, nbins = 150, nolost=1, ref=1)

    if configuration == 1:  # Be 200 um
        added_reflectivity = 1.0
        bl = set_beamline_twoCRLs(configuration=1)
    if configuration == 2:  # Be 200 um
        added_reflectivity = 1.0
        bl = set_beamline_twoCRLs(configuration=2)
    if configuration == 3:  # Be 200 um
        added_reflectivity = 1.0
        bl = set_beamline_twoCRLs(configuration=3)
    if configuration == 4: 
        added_reflectivity = 0.7
        bl = set_beamline_ml(configuration=4)
    if configuration == 5: 
        added_reflectivity = 0.7
        bl = set_beamline_ml(configuration=5)
    if configuration == 6: 
        added_reflectivity = 0.5
        bl = set_beamline_kb(configuration=6)
    if configuration == 7: 
        added_reflectivity = 0.5
        bl = set_beamline_kb(configuration=7)


    #trace
    bl.dump_systemfile()
    beam.traceCompoundOE(bl,write_start_files=1,write_end_files=1,\
        write_star_files=0,write_mirr_files=0)

    print(bl.info())
    beam.write("final.01")
    

    tk_h = beam.histo1(col=1, nbins = 150, nolost=1, ref=1)
    tk_v = beam.histo1(col=3, nbins = 150, nolost=1, ref=1)
    tk_e = beam.histo1(col=11, nbins = 150, nolost=1, ref=1)
 
    txt = ""
    #get image sizes:
    #txt += "\n\n =================================\n"
    #txt += "Ring section: %s, undulator: %s, bl config: %d\n\n"%(ring,undulator,configuration)
    txt += "==Configuration:%d, Ring section:%s ==\n\n"%(configuration,ring)
    txt += "| Source size (from setup) | H: %6.3f um | V: %6.3f um  \n"%(\
        1e4*2.35*src.SIGMAX, 1e4*2.35*src.SIGMAZ)
    txt += "| Source size (from histo1) | H: %6.3f um | V: %6.3f um \n"%(\
        1e4*tk_s_h['fwhm'], 1e4*tk_s_v['fwhm'])
    txt += "| Spot size (from stDev)    | H: %6.3f um | V: %6.3f um \n"%(\
        1e4*2.35*beam.get_standard_deviation(1,nolost=1,ref=1), \
        1e4*2.35*beam.get_standard_deviation(3,nolost=1,ref=1))
    txt += "| Spot size (from histo1)   | H: %6.3f um | V: %6.3f um \n"%(\
        1e4*tk_h['fwhm'], 1e4*tk_v['fwhm'])
    txt += "| Demagnification (real)    | H: %6.3f | V: %6.3f    \n"%(\
        tk_s_h['fwhm']/tk_h['fwhm'], tk_s_v['fwhm']/tk_v['fwhm'])
    if flux_1ev_bandwidth != 1.0:
        txt += "| N photons at source in 0.1%% bw | F0 = %.2g \n"%( flux_1ev_bandwidth*(photon_energy_ev*1e-3))
        txt += "| N photons at source in 1eV bw | F1 = %.2g \n"%( flux_1ev_bandwidth)

    e0 = photon_energy_bandwidth
    e1 = tk_e['fwhm']
    i0 = float(beam.nrays())
    i1 = beam.intensity(nolost=1)
    i0g = float(beam.nrays(nolost=1))

    if mono: txt += "| Energy bandwdith | at source: DE0=%.3f eV | at sample: DE1=%.3f eV\n"%(e0,e1)
    txt += "| Good rays | at source: N0=%d | at sample: N1=%d (%.3f %%) \n"%(int(i0),int(i0g),100.0*i0g/i0)
    txt += "| Intensity | at source: I0=%.0f | at sample: I1=%.3f (%.3f %%) \n"%(i0,i1,100*i1/i0)
    txt += "| Transmittivity excluding loses | I1/N1= | %.3f (%.3f %%) \n"%(i1/i0g,100*i1/i0g)
    if mono: txt += "| Transmittivity per eV | T1=(I1/DE1)/(I0/DE0)= | %.6f \n"%( (i1/e1)/(i0/e0) )

    if flux_1ev_bandwidth != 1.0:
        if mono: 
            if added_reflectivity == 1:
                txt += "| Integrated flux | F1 * T1 * DE1= | %.2g  \n"%( flux_1ev_bandwidth*(i1/e1)/(i0/e0)*e1 )
            else: 
                txt += "| Integrated flux | F1 * T1 * DE1 * %.1f= | %.2g  \n"%( added_reflectivity,added_reflectivity*flux_1ev_bandwidth*(i1/e1)/(i0/e0)*e1 )
    else:
        if mono: txt += "| Integrated flux | F1 * %.6f \n"%( (i1/e1)/(i0/e0)*e1 )

    return txt

if __name__ == "__main__":
    photon_energy_ev = 14200.0
    photon_energy_bandwidth = 20.0
    undulator = "U20"
    mono = 1

    #comment this for loop calculation
    #configuration = 7
    #ring = "ESRF_LB_OB"
    #ring = "ESRF_NEW_OB"

    ihit = 0
    try:
        configuration
    except NameError:
        print("Loop calculation")
        txt_all = ""
        configurations = (1,3,4,5,6,7)
        #configurations = (6,7)
        rings = ("ESRF_LB_OB","ESRF_NEW_OB")
        for configuration in configurations:
            for ring in rings:
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",ring,configuration)
                txt_all += "\n"
                txt = main()
                txt_all += txt

                tmp = Shadow.Beam()
                tmp.load("final.01")
                vcenter = numpy.mean(tmp.getshonecol(3,nolost=1))
                Shadow.ShadowTools.plotxy_gnuplot(tmp,1,3,nolost=1,ref=1,nbins=101,xrange=[-.001,.001],yrange=[vcenter-0.001,vcenter+0.001],ps=1,viewer="ls ",title="configuration:%d %s"%(configuration,ring))
                file_png = "configuration%1d%s.png"%(configuration,ring)
                os.system("convert -rotate 90 plotxy.ps %s"%(file_png))
                txt_all += "\n\n [%s] \n\n"%(file_png)
               
            txt_all += "============================\n\n"
    
        file_out = 'comparison_id23-2.t2t'
        f = open(file_out,"w")
        f.write(txt_all)
        f.close()
        print("File written to disk: ",file_out)

        # note that summary_comparison_id23-2.t2t has been written by hand
        # and include the results of comparison_id23-2.t2t automatically
        # generated
        os.system("./txt2tags -t html summary_comparison_id23-2.t2t")
        os.system("open summary_comparison_id23-2.html")

        #os.system("./txt2tags -t tex results_esrf_summary.txt")
        #os.system("pdflatex results_esrf_summary.tex")
        #os.system("okular results_esrf_summary.pdf")


    else:
        print("Single calculation, configuration=%d\n"%(configuration))
        ihit = 1
        txt = main()
        print(txt)


