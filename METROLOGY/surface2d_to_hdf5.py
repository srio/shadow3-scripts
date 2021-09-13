import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from srxraylib.plot.gol import plot
from oasys.util.oasys_util import write_surface_file
from srxraylib.metrology.profiles_simulation import slopes



# def transform_data(file_name):
#
#     """First chapuza to create a file similar to FEA"""
#
#     df = pd.read_csv(file_name, sep=';', header=None, skiprows=23)
#     # new columns #
#     df.columns = ['x(m)', 'y(m)', 'uz(m)']
#
#     new_col = ['z(m)','ux(m)','uy(m)']
#     # adding zeros for each new column
#     for col in new_col:
#         df[col] = 0.0
#
#     # reordering the columns #
#
#     cols = df.columns.tolist()
#
#     # order to be like FEA ESRF #
#     cols = cols[:2]+cols[3:4]+cols[-2:]+cols[2:3]
#
#     df = df[cols]
#
#     return df
#
# def get_line(file_name, row = 'central'):
#     """Function to get a profile file for a given Sagittal line
#     of a mirror 2D measurements"""
#
#     df = pd.read_csv(file_name, sep=';', header=None, skiprows=23)
#
#     df.columns = ['x(m)', 'y(m)', 'z(m)']
#
#     #sagittal_rows = df[df.duplicated(['y(m)'])]
#     #print(sagittal_rows)
#
#     rows_shape = df.pivot_table(columns=['y(m)'], aggfunc='size')
#
#     n_rows = rows_shape.size
#
#     if row == 'central':
#         n = int(n_rows/2)
#     elif (isinstance(row, int) == True) and (row < n_rows):
#         n = row
#     else:
#         raise RuntimeError(f'ERROR: {row} is not an integer number or is higher than the number of rows {n_rows}')
#
#     #print(rows_shape.index[n])
#
#     sub_df = df[df['y(m)'] == rows_shape.index[n]]
#
#     return sub_df
    
def get_shadow_h5(file_name):
    """Function to get an h5 file with OASYS structure
    from 2D measurements """
    
    df = pd.read_csv(file_name, sep=';', header=None, comment='#', skiprows=1)
    
    df.columns = ['x(m)', 'y(m)', 'z(m)']

    # this part is to get the ordinates and the number of abscissas for each
    rows_shape = df.pivot_table(columns=['y(m)'], aggfunc='size')

    #print(rows_shape)

    #n_rows = rows_shape.size
    
    #print(n_rows)
    
    x_coors = []
    x_mins = []
    x_maxs = []
    z_heights = []
    
    for i,y in enumerate(rows_shape.index):
        sub_df = df[df['y(m)'] == y]
        x_coors.append(np.array(sub_df['x(m)']))
        x_mins.append(x_coors[i][0])
        x_maxs.append(x_coors[i][-1])
        z_heights.append(np.array(sub_df['z(m)']))

    # checking that all coordinates along the mirror have the same steps #
    if (all(x==x_mins[0] for x in x_mins)) and (all(x==x_maxs[0] for x in x_maxs)):
        print("All elements in x_coors are the same")
        x = x_coors[0]
        y = rows_shape.index
    else:
        #TODO: define coordinates along the mirror and interpolate all#
        #z for all y coord #
        pass 
    
    #print(z_heights)
        
    return np.array(x), np.array(y), np.array(z_heights)
    
# def app_gaussian(z, sigma_0= 10, sigma_1 = 10):
#
#    """Copy paste of Manolos filtering function"""
#
#    filtered_z = gaussian_filter(z, (sigma_0,sigma_1), order=0, output=None, mode='nearest', cval=0.0, truncate=4.0)
#
#    return filtered_z
#
# def scale_profile(surface, factor):
#     """Brief function just to rescale the full surface"""
#     z2 = np.copy(surface)
#     z2 *= factor
#
#     return z2
#
#
# def detrend_best_circle(x,y,z,fitting_domain_ratio=0.5, plotting = False):
#
#     """Almost copy paste of Manolos detrend best circle function"""
#
#     xm = x.copy()
#     zm = z[y.size//2,:]
#     print(f'Medium line at {y.size//2}')
#     zm.shape = -1
#
#     icut = np.argwhere(np.abs(xm) <= fitting_domain_ratio)
#     if len(icut) <=5:
#         raise Exception("Not enough points for fitting.")
#
#     xcut = xm[icut]
#     #print(len(xm),len(xcut))
#     zmcut = zm[icut]
#
#     #print(len(zm), len(zmcut))
#
#     xcut.shape = -1
#     zmcut.shape = -1
#
#     if plotting:
#         plot(xm, zm, legend=["original"])
#
#     print( np.argwhere(np.isnan(z)))
#     print("Fitting interval: [%g,%g] (using %d points)" % (xcut[0],xcut[-1],xcut.size))
#
#     coeff = np.polyfit(xcut, np.gradient(zmcut,xcut), deg=1)
#
#     # # zfit = coeff[0] * xm  + coeff[1]
#     radius = 1 / coeff[0]
#     #print("Detrending straight line on sloped (axis=%d): zfit = %g * coordinate + %g " % (axis, coeff[1], coeff[0]))
#     print("Radius of curvature: %g m" % (1.0 / coeff[0]))
#
#     if radius >= 0:
#         zfit = radius - np.sqrt(radius ** 2 - xm ** 2)
#     else:
#         zfit = radius + np.sqrt(radius ** 2 - xm ** 2)
#     if plotting:
#         plot(xm, zfit, legend=["fit"])
#     #plot(xcut, zmcut, xm, zfit, legend=["cut","fit"])
#
#     #print(len(zfit))
#
#         plot(xm, zm-zfit, legend=["detrended"])
#
#     for i in range(z.shape[0]):
#         z[i,:] -= zfit
#
#
#     nx, ny = z.shape
#     z = z - (z[nx//2,ny//2])
#
#     # print(f" Slope error is {round(z[:, 0].std(), 6)}")
#
#     return xm, z
    
def plot2d(x,y,data):
		
    plt.pcolormesh(x,y,data, cmap=plt.cm.viridis)
	
    plt.colorbar().ax.tick_params(axis='y',labelsize=12)

    plt.ylabel("Vertical [mm]",fontsize=12)
    plt.xlabel("Horizontal [mm]",fontsize=12)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

if __name__ == '__main__':


    file_name = 'ring256_TypbeB_F127001_frontside_ontissue_meas2__avg_2D.txt'
    x, y, z = get_shadow_h5(file_name)
    print(z.shape, x.shape, y.shape, z.min(), z.max())

    from srxraylib.plot.gol import plot_image
    plot_image(z*1e6, y*1e3, x*1e3, aspect="auto")
    
    # x,z  = detrend_best_circle(x,y,z,fitting_domain_ratio=0.5, plotting=True)
    #
    # print(z.shape)
    # #plot2d(x,y,z)
    #
    # z2 = app_gaussian(z, sigma_0= 6, sigma_1 = 2)
    #
    # z3 =  scale_profile(z2,1)
    #
    # #plot2d(x,y,z)
    slp = slopes(z, y, x, silent=0, return_only_rms=0)
    #
    # slp_y = np.round(slp[1][1]*1e6, 3)

    output_filename = f'ring256.h5'

    # plot(x,z[y.size//2,:],x,z[y.size//2,:],legend=["detrended","Gauss_filtered"])
    #
    # plot(x,np.gradient(z[y.size//2,:],x), legend=["Slope errors"])

    write_surface_file(z.T, y, x, output_filename, overwrite=True)
    print("write_h5_surface: File for OASYS " + output_filename + " written to disk.")

    print(">>>>>", z.T.shape, y.shape, x.shape,)