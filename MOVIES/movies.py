from orangecontrib.comsyl.util.CompactAFReader import CompactAFReader
import numpy as np
import matplotlib.pylab as plt


def plot_modes(ring='ebs',test=False):

    #
    # mogrify -resize 600x400  *.png
    # convert *.png -delay 10 -morph 10 %05d.png
    # ffmpeg -r 25  -i %05d.png output.mp4





    if ring == 'ebs':
        filename="/scisoft/data/srio/COMSYL/ID16/id16s_ebs_u18_1400mm_1h_new_s1.0.npz"
        max_modes = 1000
        aspect = None
    elif ring == 'hb':
        filename="/scisoft/data/srio/COMSYL/ID16/id16s_hb_u18_1400mm_1h_s1.0.h5"
        max_modes = 4000
        aspect = 'auto'
    elif ring == 'ebs_id09':
        file_dir = "/scisoft/users/srio/COMSYLD/comsyl/comsyl/calculations/"
        filename=file_dir+"id09_ebs_u17_2000mm_1h_s1.5.npy"
        max_modes = 800
        aspect = None
    elif ring == 'lb_id09':
        file_dir = "/scisoft/users/srio/COMSYLD/comsyl/comsyl/calculations/"
        filename=file_dir+"id09_lowbeta_u17_2000mm_1h_s2.5.npy"
        max_modes = 1500
        aspect = 'auto'
    else:
        raise Exception('Not implemented')

    af  = CompactAFReader.initialize_from_file(filename)

    x = 1e6 * af.x_coordinates()
    y = 1e6 * af.y_coordinates()

    nx = 12
    nfactor = 2./3

    max_mode_index = 24




    modes_to_plot = np.concatenate((np.arange(0,max_mode_index,1),np.arange(100,max_modes,50)))
    print (">>>>",modes_to_plot)

    if test:
        modes_to_plot = [0]

    for i in modes_to_plot:

        fig = plt.figure(figsize=(nx,int(nx*nfactor)))

        if ring == 'ebs':
            left, width    = 0.1, 0.3
            bottom, height = 0.2, 0.3
            rect_scatter = [left, bottom, width, height]
            axScatter = plt.axes(rect_scatter)
            plt.title("Cumulated occupation")


            left, width    = 0.1,    0.4
            bottom, height = 0.65,   0.4*nfactor
            rect_scatter = [left, bottom, width, height]
            ayScatter = plt.axes(rect_scatter)
            plt.title("Single mode")


            left, width    = 0.475, 0.5
            bottom, height = 0.2, 0.5*nfactor
            rect_scatter = [left, bottom, width, height]
            azScatter = plt.axes(rect_scatter)
            plt.title("Cumulated Spectral Density")


            left, width    = 0.575, 0.4
            bottom, height = 0.65,  0.4*nfactor
            rect_scatter = [left, bottom, width, height]
            akScatter = plt.axes(rect_scatter)
            plt.title("Spectral Density")
        elif ring == 'ebs_id09':
            left, width    = 0.1, 0.3
            bottom, height = 0.2, 0.3
            rect_scatter = [left, bottom, width, height]
            axScatter = plt.axes(rect_scatter)
            plt.title("Cumulated occupation")


            left, width    = 0.1,    0.4
            bottom, height = 0.65,   0.4*nfactor
            rect_scatter = [left, bottom, width, height]
            ayScatter = plt.axes(rect_scatter)
            plt.title("Single mode")


            left, width    = 0.475, 0.5
            bottom, height = 0.2, 0.5*nfactor
            rect_scatter = [left, bottom, width, height]
            azScatter = plt.axes(rect_scatter)
            plt.title("Cumulated Spectral Density")


            left, width    = 0.575, 0.4
            bottom, height = 0.65,  0.4*nfactor
            rect_scatter = [left, bottom, width, height]
            akScatter = plt.axes(rect_scatter)
            plt.title("Spectral Density")
        elif ring == 'hb':

            left, width    = 0.1, 0.3
            bottom, height = 0.2, 0.3
            rect_scatter = [left, bottom, width, height]
            axScatter = plt.axes(rect_scatter)
            plt.title("Cumulated occupation")



            left, width    = 0.1,    0.4
            bottom, height = 0.65,   0.4*nfactor
            rect_scatter = [left, bottom, width, height]
            ayScatter = plt.axes(rect_scatter)
            plt.title("Single mode")


            left, width    = 0.475, 0.5
            bottom, height = 0.15, 0.5*nfactor
            rect_scatter = [left, bottom, width, height]
            azScatter = plt.axes(rect_scatter)
            plt.title("Cumulated Spectral Density")

            left, width    = 0.575, 0.4
            bottom, height = 0.65,  0.4*nfactor
            rect_scatter = [left, bottom, width, height]
            akScatter = plt.axes(rect_scatter)
            plt.title("Spectral Density")

        else:
            raise Exception('Not implemented')





        x_ebs = np.arange(af.number_modes())
        y_ebs = np.cumsum(np.abs(af.occupation_array()))


        # cumulated spectrum
        axScatter.plot(x_ebs,y_ebs)
        axScatter.plot(x_ebs[i],y_ebs[i],"o")
        axScatter.set_ylabel("Cumulated occupation")
        axScatter.set_xlabel("Mode index")


        # mode i
        ayScatter.imshow(np.abs(af.mode(i)).T,origin='lower',extent=[x[0],x[-1],y[0],y[-1]],aspect=aspect)
        ayScatter.set_ylabel("X [$\mu$m]")
        ayScatter.set_xlabel("Y [$\mu$m]")


        # CSD
        azScatter.imshow(af.intensity_from_modes(max_mode_index=i).T,origin='lower',extent=[x[0],x[-1],y[0],y[-1]],aspect=aspect)
        azScatter.set_xlabel("X [$\mu$m]")
        azScatter.set_ylabel("Y [$\mu$m]")

        # Tital CSD
        akScatter.imshow(np.abs(af.spectral_density()).T,origin='lower',extent=[x[0],x[-1],y[0],y[-1]],aspect=aspect)
        akScatter.set_xlabel("X [$\mu$m]")
        akScatter.set_ylabel("Y [$\mu$m]")


        if not test:
            if ring == 'ebs':
                fileout = "MOVIE_EBS/BACKUP/movie%04d.png"%i
            elif ring == 'hb':
                fileout = "MOVIE_HB/BACKUP/movie%04d.png"%i
            elif ring == 'ebs_id09':
                fileout = "MOVIE_EBS/BACKUP/movie%04d.png"%i
            elif ring == 'lb_id09':
                fileout = "MOVIE_LB/BACKUP/movie%04d.png"%i
            else:
                raise Exception('Not implemented')
            plt.savefig(fileout)
            print("File written to disk: %s"%fileout)
            plt.close('all')

    if test:
        plt.show()




if __name__ == "__main__":

    plot_modes(ring='ebs_id09')

    #plot_modes(ring='hb',test=True)






