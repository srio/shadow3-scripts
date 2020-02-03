from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy as np

from srxraylib.plot.gol import plot_image, plot
import matplotlib.pylab as plt


def up_to_mode(files,titles,thresholds=[0.95,0.75,0.50]):


    fup = open("up_to_mode_number.dat",'w')
    for i,file in enumerate(files):
        af  = CompactAFReader.initialize_from_file(file)
        x_ebs = np.arange(af.number_modes())
        y_ebs = np.cumsum(np.abs(af.occupation_array()))



        print("\n>>>>>>>>>>>>>>>>: ",titles[i])
        for threshold in thresholds:
            # threshold = 0.95O
            i_good = np.where(y_ebs > threshold)
            print(">>>>>>>>>>>>>>>> mode up to %d percent: %d"%(threshold*100, x_ebs[i_good[0][0]]))

        plt.plot(x_ebs,y_ebs,label=titles[i])

        filename = "cumulated_occupation_%s.dat"%(titles[i])
        f = open(filename,'w')
        for j,x in enumerate(x_ebs):
            f.write("%g  %g\n"%(x,y_ebs[j]))
        f.close()
        print("File written to disk: %s"%filename)

    fup.close()

    ax = plt.subplot(111)
    ax.legend(bbox_to_anchor=None)

    plt.ylabel("Cumulated occupation")
    plt.xlabel("Mode index")
    plt.savefig("cumulated_occupation.eps")
    plt.show()


def find_ebs_file_among_all_u18_2m():

    path = "/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/"

    list = [
        "cd_01_u18_2m_s1.4.npy",
        "cd_01_u18_2m_s2.0.npy",
        "cd_01_u18_2m_s2.5.npy",
        "cd_01_u18_2m_s3.0.npy",
        "cd_01_u18_2m_s3.5.npy",
        "cd_05_u18_2m_s2.0.npy",
        "cd_05_u18_2m_s2.5.npy",
        "cd_05_u18_2m_s3.0.npy",
        "cd_05_u18_2m_s3.5.npy",
        "cd_12_u18_2m_s2.0.npy",
        "cd_12_u18_2m_s2.5.npy",
        "cd_12_u18_2m_s3.0.npy",
        "cd_12_u18_2m_s3.5.npy",
        "cd_16_u18_2m_s2.0.npy",
        "cd_16_u18_2m_s2.5.npy",
        "cd_16_u18_2m_s3.0.npy",
        "cd_16_u18_2m_s3.5.npy",
        "cd_1_u18_2m_s2.0.npy",
        "cd_1_u18_2m_s2.5.npy",
        "cd_1_u18_2m_s3.0.npy",
        "cd_1_u18_2m_s3.5.npy",
        "cd_20_u18_2m_s2.0.npy",
        "cd_20_u18_2m_s2.5.npy",
        "cd_20_u18_2m_s3.0.npy",
        "cd_20_u18_2m_s3.5.npy",
        "cd_24_u18_2m_s2.0.npy",
        "cd_24_u18_2m_s2.5.npy",
        "cd_24_u18_2m_s3.0.npy",
        "cd_24_u18_2m_s3.5.npy",
        "cd_28_u18_2m_s2.0.npy",
        "cd_28_u18_2m_s2.5.npy",
        "cd_28_u18_2m_s3.0.npy",
        "cd_28_u18_2m_s3.5.npy",
        "cd_4_u18_2m_s2.0.npy",
        "cd_4_u18_2m_s2.5.npy",
        "cd_4_u18_2m_s3.0.npy",
        "cd_4_u18_2m_s3.5.npy",
        "cd_8_u18_2m_s2.0.npy",
        "cd_8_u18_2m_s2.5.npy",
        "cd_8_u18_2m_s3.0.npy",
        "cd_8_u18_2m_s3.5.npy",
        "cdt_01_u18_2m_s2.0.npy",
        "cdt_01_u18_2m_s2.5.npy",
        "cdt_01_u18_2m_s3.0.npy",
        "cdt_01_u18_2m_s3.5.npy",
        "cdt_05_u18_2m_s2.0.npy",
        "cdt_05_u18_2m_s2.5.npy",
        "cdt_05_u18_2m_s3.0.npy",
        "cdt_05_u18_2m_s3.5.npy",
        "cdt_12_u18_2m_s2.0.npy",
        "cdt_12_u18_2m_s2.5.npy",
        "cdt_12_u18_2m_s3.0.npy",
        "cdt_12_u18_2m_s3.5.npy",
        "cdt_16_u18_2m_s2.0.npy",
        "cdt_16_u18_2m_s2.5.npy",
        "cdt_16_u18_2m_s3.0.npy",
        "cdt_16_u18_2m_s3.5.npy",
        "cdt_1_u18_2m_s2.0.npy",
        "cdt_1_u18_2m_s2.5.npy",
        "cdt_1_u18_2m_s3.0.npy",
        "cdt_1_u18_2m_s3.5.npy",
        "cdt_20_u18_2m_s2.0.npy",
        "cdt_20_u18_2m_s2.5.npy",
        "cdt_20_u18_2m_s3.0.npy",
        "cdt_20_u18_2m_s3.5.npy",
        "cdt_24_u18_2m_s2.0.npy",
        "cdt_24_u18_2m_s2.5.npy",
        "cdt_24_u18_2m_s3.0.npy",
        "cdt_24_u18_2m_s3.5.npy",
        "cdt_28_u18_2m_s2.0.npy",
        "cdt_28_u18_2m_s2.5.npy",
        "cdt_28_u18_2m_s3.0.npy",
        "cdt_28_u18_2m_s3.5.npy",
        "cdt_4_u18_2m_s2.0.npy",
        "cdt_4_u18_2m_s2.5.npy",
        "cdt_4_u18_2m_s3.0.npy",
        "cdt_4_u18_2m_s3.5.npy",
        "cdt_8_u18_2m_s2.0.npy",
        "cdt_8_u18_2m_s2.5.npy",
        "cdt_8_u18_2m_s3.0.npy",
        "cdt_8_u18_2m_s3.5.npy",
        "cl_high_beta_u18_2m_1h_s1.4.npy",
        "cl_high_beta_u18_2m_1h_s1.5.npy",
        "cl_high_beta_u18_2m_1h_s1.6.npy",
        "cl_high_beta_u18_2m_1h_s1.7.npy",
        "cl_high_beta_u18_2m_1h_s1.8.npy",
        "cl_high_beta_u18_2m_1h_s1.9.npy",
        "cl_high_beta_u18_2m_1h_s2.0.npy",
        "cl_low_beta_u18_2m_1h_s6.5.npy",
        "cm_delta_u18_2m_1h_matrix_s0.5.npy",
        "cm_delta_u18_2m_1h_matrix_s1.0.npy",
        "cm_delta_u18_2m_1h_matrix_s1.5.npy",
        "cm_delta_u18_2m_1h_matrix_s2.0.npy",
        "cm_delta_u18_2m_1h_tsnc_s1.0.npy",
        "cm_delta_u18_2m_1h_tsnc_s2.0.npy",
        "cm_delta_u18_2m_1h_tsnc_s3.0.npy",
        "cm_delta_u18_2m_1h_tsnc_s4.0.npy",
        "cm_delta_u18_2m_1h_ts_s1.0.npy",
        "cm_delta_u18_2m_1h_ts_s2.0.npy",
        "cm_delta_u18_2m_1h_ts_s3.0.npy",
        "cm_delta_u18_2m_1h_ts_s4.0.npy",
        "cm_delta_u18_2m_1h_ts_s5.0.npy",
        "cm_delta_u18_2m_1h_ts_s6.0.npy",
        "cm_delta_u18_2m_1h_ts_s7.0.npy",
        "cm_delta_u18_2m_1h_ts_s8.0.npy",
        "cm_new_alpha_u18_2m_1h_matrix_s1.5.npy",
        "cm_new_alpha_u18_2m_1h_tsnc_s1.5.npy",
        "cm_new_alpha_u18_2m_1h_tsnc_s2.5.npy",
        "cm_new_res_u18_2m_1h_matrix_s2.5.npy",
        "cm_new_res_u18_2m_1h_tsnc_s2.5.npy",
        "cm_new_res_u18_2m_1h_ts_s2.5.npy",
        "cm_new_u18_2m_1h_matrix_s2.5.npy",
        "cm_new_u18_2m_1h_ts_s2.5.npy",
        "cp_new_u18_2m_1h_ts_s2.5.npy",
        "cs_new_alpha_u18_2m_1h_s1.8.npy",
        "cs_new_u18_2m_1h_s2.5.npy",
        "cu_new_u18_2m_1h_s1.0.npy",
        "cu_new_u18_2m_1h_s1.5.npy",
        "cu_new_u18_2m_1h_s2.0.npy",
        "cu_new_u18_2m_1h_s2.5.npy",
        "cu_new_u18_2m_1h_s3.0.npy",
        "cu_new_u18_2m_1h_s3.5.npy",
        "cu_new_u18_2m_1h_s4.0.npy",
        "cu_new_u18_2m_3h_s1.0.npy",
        "cu_new_u18_2m_3h_s1.5.npy",
        "cu_new_u18_2m_3h_s2.0.npy",
        "cu_new_u18_2m_3h_s2.5.npy",
        "cu_new_u18_2m_3h_s3.0.npy",
        "cu_new_u18_2m_3h_s3.5.npy",
        "cu_new_u18_2m_3h_s4.0.npy",
        "cu_new_u18_2m_5h_s0.2.npy",
        "cu_new_u18_2m_5h_s0.4.npy",
        "cu_new_u18_2m_5h_s0.6.npy",
        "cu_new_u18_2m_5h_s0.8.npy",
        "cu_new_u18_2m_5h_s1.0.npy",
        "cu_new_u18_2m_5h_s1.2.npy",
        "cu_new_u18_2m_5h_s1.4.npy",
        "cu_new_u18_2m_5h_s1.6.npy",
        "delta_u18_2m_1h_s1.0.npy",
        "delta_u18_2m_1h_s2.0.npy",
        "delta_u18_2m_1h_s3.0.npy",
        "delta_u18_2m_1h_s4.0.npy",
        "delta_u18_2m_1h_s5.0.npy",
        "gel_01_u18_2m_s2.5.npy",
        "gel_05_u18_2m_s2.5.npy",
        "gel_12_u18_2m_s2.5.npy",
        "gel_16_u18_2m_s2.5.npy",
        "gel_1_u18_2m_s2.5.npy",
        "gel_20_u18_2m_s2.5.npy",
        "gel_4_u18_2m_s2.5.npy",
        "gel_8_u18_2m_s2.5.npy",
        "gel_big2_u18_2m_s1.0.npy",
        "gel_big_u18_2m_s1.0.npy",
        "hf_new_u18_2m_1h_s0.5.npy",
        "hf_new_u18_2m_1h_s1.0.npy",
        "hf_new_u18_2m_1h_s1.5.npy",
        "hf_new_u18_2m_1h_s2.0.npy",
        "hf_new_u18_2m_1h_s2.5.npy",
        "hf_new_u18_2m_1h_s3.0.npy",
        "hf_new_u18_2m_1h_s3.5.npy",
        "hf_new_u18_2m_1h_s4.0.npy",
        "hf_new_u18_2m_3h_s1.0.npy",
        "hf_new_u18_2m_3h_s1.5.npy",
        "hf_new_u18_2m_3h_s2.0.npy",
        "hf_new_u18_2m_3h_s2.5.npy",
        "hf_new_u18_2m_3h_s3.0.npy",
        "hf_new_u18_2m_3h_s3.5.npy",
        "hf_new_u18_2m_3h_s4.0.npy",
        "low_beta_res_u18_2m_1h_ts_s4.5.npy",
        "low_beta_res_u18_2m_1h_ts_s5.0.npy",
        "low_beta_res_u18_2m_1h_ts_s5.1.npy",
        "low_beta_res_u18_2m_1h_ts_s6.0.npy",
        "low_beta_res_u18_2m_1h_ts_s6.5.npy",
        "low_beta_res_u18_2m_1h_ts_s8.0.npy",
        "new_res_u18_2m_1h_matrix_s2.0.npy",
        "new_res_u18_2m_1h_matrix_s2.5.npy",
        "new_res_u18_2m_1h_tsnc_s2.0.npy",
        "new_res_u18_2m_1h_ts_s2.0.npy",
        "new_res_u18_2m_1h_ts_s2.5.npy",
        "new_u18_2m_1h_matrix_s2.0.npy",
        "new_u18_2m_1h_s0.1.npy",
        "new_u18_2m_1h_s1.0.npy",
        "new_u18_2m_1h_ts_s2.0.npy",
        "ns_cl_high_beta_u18_2m_1h_s1.6.npy",
        "test_cm_new_u18_2m_1h_ts_s2.5.npy",]

    for file in list:
        print("%s"%file)
        af  = CompactAFReader.initialize_from_file(path+"/"+file)
        y_ebs = np.cumsum(np.abs(af.occupation_array()))
        if np.abs(y_ebs[0]-7.015756401120343266e-02) < 1e-4:
            print(">>>>>>>>>>>>>>>>>>found: %s"%file)





if __name__ == "__main__":

    # find_ebs_file_among_all_u18_2m()



    # filename_ebs="/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cs_new_u18_2m_1h_s2.5.npy"
    # filename_hb="/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_high_beta_u18_2m_1h_s1.4.npy"
    # filename_lb="/scisoft/users/glass/Documents/sources/Orange-SRW/comsyl/calculations/cl_low_beta_u18_2m_1h_s6.5.npy"
    # files = [filename_ebs,filename_lb,filename_hb]
    # titles=["ebs","low_beta","high_beta"]
    # up_to_mode(files,titles)





    filename_ebs="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npy"
    filename_hb ="/scisoft/data/srio/COMSYL/ID16/id16s_hb_u18_1400mm_1h_s1.0.npy"
    files = [filename_ebs,filename_hb]
    titles=["id16a_ebs","id16a_high_beta"]
    up_to_mode(files,titles,thresholds=[0.85,0.75,0.50])


