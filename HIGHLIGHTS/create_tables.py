#
# runs using modified wofry propagator and adapted for numpy wavefronts ans pickle
#



import numpy


# from srxraylib.plot.gol import plot_image, plot
# from wofry.propagator.propagator import PropagationElements, PropagationParameters
# from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D






if __name__ == "__main__":

    from comsyl.waveoptics.WofrySuperBeamline import WofrySuperBeamline
    from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader



    filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    # filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy" # OK LB
    # filename = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s2.0.npy"

    #
    # load CSD
    #

    af  = CompactAFReader.initialize_from_file(filename)
    intensity_full = af.total_intensity()


    # x = af.x_coordinates()
    # y = af.y_coordinates()
    # cumulated_occupation = af.cumulated_occupation_array()
    # occupation = af.occupation_array()



    slitV = [30,25,20,15,10,5, 5] # [25,20,15,10,5] #
    slitH = [25,25,20,15,10,10,5] # [25,20,15,10,5] #



    #
    # intensities / spectral density fraction
    #
    if True:
        f = open("propagation.txt",'w')
        f.close()

        for j in range(len(slitV)):

            file_propagated = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/NICE/propagation_wofry_%dx%d/beamline_cs_new_u18_2m_1h_s2.5.npz"%(slitV[j],slitH[j])

            af_prop  = CompactAFReader.initialize_from_file(file_propagated)
            intensity_prop = af_prop.total_intensity()
            intensity_prop = af_prop.total_intensity_from_modes()
            ratio = intensity_prop/intensity_full

            print("slit H:%d, V:%d: intensity(full): %g intensity(cut): %g ratio: %g "%(
                slitH[j],slitV[j],intensity_full,intensity_prop,ratio,))

            # cum0 = 0.0
            # cum1 = 0.0
            # for i in range(af_prop.number_modes()):
            #     # print(af.occupation(i),af_prop.occupation(i))
            #     # print(af.mode_intensity(i),af_prop.mode_intensity(i))
            #     cum0 += numpy.abs(af.eigenvalue(i))
            #     cum1 += numpy.abs(af_prop.eigenvalue(i) * af_prop.mode_intensity(i))
            #     print(i,af.eigenvalue(i),af_prop.mode_intensity(i),cum0,cum1)
            # print("FINAL: ",cum0,cum1,cum1/cum0)

            # print("\n\n\n\n\n\n")



            f = open("propagation.txt",'a')
            f.write('"%d x %d" %f \n'%(slitV[j],slitH[j],ratio))
            f.close
        print("File written to disk: propagation.txt")


    #
    # coherent fraction and occupation after rediagonalization
    #

    if True:
        f = open("rediagonalized.txt",'w')
        f.close()

        OCCUPATION = []
        title = "mode "
        for j in range(len(slitV)):

            file_propagated = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/NICE/rediagonalized_%sx%d.npz"%(slitV[j],slitH[j])
            print("Finding file: ",file_propagated)

            af_prop  = CompactAFReader.initialize_from_file(file_propagated)

            intensity_prop = af_prop.total_intensity()
            intensity_prop = af_prop.total_intensity_from_modes()
            ratio = intensity_prop/intensity_full

            print("slit H:%d, V:%d: intensity(full): %g intensity(cut): %g ratio: %g "%(
                slitH[j],slitV[j],intensity_full,intensity_prop,ratio,))

            cum0 = 0.0
            cum1 = 0.0

            print(af_prop.occupation_array())
            OCCUPATION.append(af_prop.occupation_array())

            f = open("rediagonalized.txt",'a')
            f.write('"%d x %d" %f \n'%(slitV[j],slitH[j],ratio))
            f.close

            title += "%sx%d    "%(slitV[j],slitH[j])


        f = open("rediagonalized.txt",'a')
        f.write("\n\n\n\n")
        f.write(title+" \n")
        for i in range(10):
            f.write("%d    "%(i))
            for j in range(len(OCCUPATION)):
                f.write("%f "%(OCCUPATION[j][i]))
            f.write("\n")

        f.close
        print("File written to disk: rediagonalized.txt ")
