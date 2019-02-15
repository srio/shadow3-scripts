#
# runs using modified wofry propagator and adapted for numpy wavefronts ans pickle
#



import numpy


# from srxraylib.plot.gol import plot_image, plot
# from wofry.propagator.propagator import PropagationElements, PropagationParameters
# from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D






if __name__ == "__main__":

    # from comsyl.waveoptics.WOFRYAdapter import CWBeamline
    from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader



    filename_ebs = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npz" # OK EBS
    filename_lb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy" # OK LB
    filename_hb = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s2.0.npy"





    rings = ["EBS","EBS","EBS","LB","LB","LB","HB","HB","HB",]
    cases = [
        "srw_EBS_25x25","wofry_EBS_25x25","wofry_EBS_15x15",
        "srw_LB_25x25","wofry_LB_25x25","wofry_LB_15x15",
        "srw_HB_25x25","wofry_HB_25x25","wofry_HB_15x15",

    ]




    #
    # coherent fraction and occupation after rediagonalization
    #

    if True:
        f = open("rediagonalized.txt",'w')
        f.close()

        OCCUPATION = []
        title = "mode "







        for i in range(len(cases)):

            #
            # load CSD
            #
            if rings[i] == "EBS":
                af  = CompactAFReader.initialize_from_file(filename_ebs)
            elif rings[i] == "HB":
                af  = CompactAFReader.initialize_from_file(filename_hb)
            elif rings[i] == "LB":
                af  = CompactAFReader.initialize_from_file(filename_lb)


            intensity_full = af.total_intensity()

            file_propagated = "/scisoft/xop2.4/extensions/shadowvui/shadow3-scripts/HIGHLIGHTS/NICE-NEW/rediagonalized_%s.npz"%cases[i]

            print("Finding file: ",file_propagated)

            af_prop  = CompactAFReader.initialize_from_file(file_propagated)

            intensity_prop = af_prop.total_intensity()
            intensity_prop = af_prop.total_intensity_from_modes()
            ratio = intensity_prop/intensity_full

            print("case: %s,  intensity(full): %g intensity(cut): %g ratio: %g "%(
                cases[i],intensity_full,intensity_prop,ratio,))

            cum0 = 0.0
            cum1 = 0.0

            print(af_prop.occupation_array())
            OCCUPATION.append(af_prop.occupation_array())

            f = open("rediagonalized.txt",'a')
            f.write("%s %f \n"%(cases[i],ratio))
            f.close

            title += "%s    "%(cases[i])


        f = open("rediagonalized.txt",'a')
        f.write("\n\n\n\n")
        f.write(title+" \n")
        for i in range(10):
            f.write("%d    "%(i))
            for j in range(len(OCCUPATION)):
                f.write("%f "%((OCCUPATION[j][i])).real)
            f.write("\n")

        f.close
        print("File written to disk: rediagonalized.txt ")
