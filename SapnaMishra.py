import traceback
import numpy as np

def write_shadow_surface(s,xx,yy,filename='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,filename=filename)
      INPUTS:
           z - 2D array of heights z(x,y)
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           filename - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """

    try:
       fs = open(filename, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+filename)
       return
    else:
        # dimensions
        fs.write( "%d  %d \n"%(xx.size,yy.size))
        # y array
        for i in range(yy.size):
            fs.write("%g  "%(yy[i]))
        fs.write("\n")
        # for each x element, the x value followed by the corresponding z(y) profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps += "%g  "%(s[i,j])
            fs.write("%g    %s \n"%(xx[i],tmps))
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+filename+" written to disk.")


def write_shadow_surface1(s,xx,yy,filename='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,filename=filename)
      INPUTS:
           z - 2D array of heights z(x,y)
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           filename - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """

    try:
       fs = open(filename, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+filename)
       return
    else:
        # dimensions
        fs.write( "%d  %d \n"%(xx.size,yy.size))
        # y array
        for i in range(yy.size):
            fs.write("%g  "%(yy[i]))
        fs.write("\n")
        # for each x element, the x value followed by the corresponding z(y) profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps += "%g  "%(s[i,j])
            fs.write("%g    %s \n"%(xx[i],tmps))
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+filename+" written to disk.")


import numpy as np

import matplotlib.pyplot as plt

import matplotlib.ticker as ticker

import pandas as pd

import sympy as sp

import math


def math_coord(B_low, B_high, n, n_p):
    B = np.linspace(B_low, B_high, n_p)  # Bragg Angle

    B_rad = np.radians(B)

    xc = ((1 / (n * np.tan(B_rad)))) ** (1 / (n - 1))

    yc = xc ** n

    r = ((1 + (n * ((xc ** (n - 1)))) ** 2) ** (3 / 2)) / (n * (n - 1) * (xc ** (n - 2)))

    r_sinb = r * np.sin(B_rad)

    xJ = xc - r_sinb * np.sin(2 * B_rad)

    yJ = yc - r_sinb * np.cos(2 * B_rad)

    return B, B_rad, xc, yc, r, r_sinb, xJ, yJ


def crystal_coord(y_offset, H, xc, yc, r, r_sinb, xJ, yJ):
    k = (H / (xc[0] - xc[-1]))

    print("k = ", k)

    new_xc1 = []
    new_yc1 = []
    new_xJ1 = []
    new_yJ1 = []
    for i in range(len(xc)):
        l = len(xc) / 2
        new_xc = xc[i] - xc[-1]
        new_yc = k * (yc[i] - yc[-1]) + y_offset
        new_xJ = xJ[i] - xc[-1]
        new_yJ = k * (yJ[i] - yc[-1]) + y_offset
        new_xc1.append(new_xc)
        new_yc1.append(new_yc)
        new_xJ1.append(new_xJ)
        new_yJ1.append(new_yJ)

    Xc = k * np.array(new_xc1)

    Yc = np.array(new_yc1)

    R = k * r

    R_sinb = k * r_sinb

    XJ = k * np.array(new_xJ1)

    YJ = np.array(new_yJ1)

    Xcr, Ycr = Xc, Yc

    XJr, YJr = XJ, YJ

    center = len(YJr) // 2

    Crys_x_ext = Xcr[0] - Xcr[-1]

    # print("Crystal X-Extent = ",Crys_x_ext)

    Crys_y_ext = Ycr[0] - Ycr[-1]

    # print("Crystal Y-Extent = ",Crys_y_ext)

    Det_len = YJr[-1] - YJr[0]

    # print("Detector Length = ",Det_len)

    Crys_len = math.sqrt((Xcr[-1] - Xcr[0]) ** 2 + (Ycr[-1] - Ycr[0]) ** 2)

    # print("Crystal length = ",Crys_len)

    Det_loc = math.sqrt((Xcr[-1] - XJr[0]) ** 2)

    # print("Detector Location = ",Det_loc)

    return Xcr, Ycr, R, R_sinb, XJr, YJr, Crys_y_ext, Det_loc


def curve_fitted_by_general_eq(Xc, Yc, Xc_new, plotting):
    # Fit a polynomial function to the data
    degree = 3  # By default we have taken 3 because that is used in all
    coefficients = np.polyfit(Xc, Yc, degree)

    # Create a symbolic variable for x
    x = sp.symbols('x')

    # Generate the equation using the fitted coefficients
    equation = sum(sp.sympify(coeff) * x ** (degree - i) for i, coeff in enumerate(coefficients))

    equation = sp.simplify(equation)

    a = coefficients[0]

    b = coefficients[1]

    c = coefficients[2]

    d = coefficients[
        3]  # remember d is just constant which will counter the effect of offset so it will be having exactly same values as offset
    dy_dx = sp.diff(equation, x)

    d2y_dx2 = sp.diff(dy_dx, x)

    if plotting == True:
        print("Equation = ", equation, '\n')  # use display to show good looking equation
        print("a = ", a, '\n')
        print("b = ", b, '\n')
        print("c = ", c, '\n')
        print("d = ", d, '\n')
        print("dy_dx = ", dy_dx, '\n')
        print("d2y_dx2 = ", d2y_dx2, '\n')

    Rc = ((1 + dy_dx ** 2) ** (3 / 2)) / (d2y_dx2)

    Yc_new = [equation.subs(x, val) for val in Xc_new]

    Rc_value = [Rc.subs(x, val) for val in Xc_new]

    Rc_value = np.array(Rc_value)

    #     print("Yc_new = ", Yc_new,'\n')

    #     print("Rc_value  = ",Rc_value,'\n')

    # In order to check plots that fitted is same as original_data or not
    if plotting == True:
        plt.scatter(Xc, Yc, label='Original Data', color='blue')
        plt.plot(Xc_new, Yc_new, label='Curve Fit', color='red')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Curve Fit Example')
        plt.legend()
        plt.grid(True)
        plt.show()
    else:
        pass

    return Yc_new, Rc_value, equation, a, b, c, d


# inputs

H = 60

y_offset = 0

B1_low = 25  # lower bragg angle

B1_high = 45  # higher bragg angle

n = 3  # power of x in equation of y = x^n here n = 3 so y = x^3

n_points = 21  # number of points between B1_low and B1_high

B, B_rad, xc, yc, r, r_sinb, xJ, yJ = math_coord(B1_low, B1_high, n, n_points)

print("input B1_low: ", B1_low)
print("input B1_high: ", B1_high)
print("input n: ", n)
print("input n_points: ", n_points)

print("input B1_low: ", B1_low)
print("input B1_high: ", B1_high)
print("input n: ", n)
print("input n_points: ", n_points)

from srxraylib.plot.gol import plot
plot(xc, yc)

Xcr, Ycr, R, R_sinb, XJr, YJr, Crys_y_ext, Det_loc = crystal_coord(y_offset, H, xc, yc, r, r_sinb, xJ, yJ)

plotting = True

Yc_new, Rc_value, equation, a, b, c, d = curve_fitted_by_general_eq(Xcr, Ycr, Xcr, plotting)
# print("Rc = ",Rc_value)

tx1 = Xcr / 10  # in cm

ty1 = Ycr / 10  # in cm

Crys_x_ext1 = tx1[0] - tx1[-1]

print("Crystal X - extent [cm] = ", Crys_x_ext1, '\n')

Crys_y_ext1 = ty1[0] - ty1[-1]

print("Crystal Y - extent [cm] = ", Crys_y_ext1, '\n')

Det_len1 = YJr[-1] / 10 - YJr[0] / 10

print("Detector Length [cm] = ", Det_len1, '\n')

Crys_len1 = math.sqrt((tx1[-1] - tx1[0]) ** 2 + (ty1[-1] - ty1[0]) ** 2)

print("Crystal Length [cm] = ", Crys_len1, '\n')

Det_loc1 = math.sqrt((tx1[-1] - XJr[0] / 10) ** 2)

print("Detector Location [cm] = ", Det_loc1)

p = np.linspace(-10, 10, len(tx1))
p = np.linspace(-15, 15, len(tx1))

# Surface

# original surface .dat file generated through below shown line
# step-1
tx1 = np.array(tx1)
tx1_sorted_indices = np.argsort(tx1)

# Z1 = np.outer(np.ones_like(ty1), ty1)
# # Z1 *= 0
# # tx1 = tx1 - tx1[tx1.size//2]
# # tx1 *= 2
# write_shadow_surface(Z1, p, tx1, "/users/srio/Oasys/cubic_surface_step1.dat")
# plot(tx1, Z1[Z1.shape[0]//2, :], title="central profile")

Z1 = np.outer(np.ones_like(ty1), ty1[tx1_sorted_indices])
tx11 = tx1[tx1_sorted_indices]
write_shadow_surface(Z1, p, tx11, "/users/srio/Oasys/cubic_surface_step1.dat")
plot(tx11, Z1[Z1.shape[0]//2, :], title="central profile")

# write file for h5
from oasys.util.oasys_util import write_surface_file
write_surface_file(Z1.T*1e-2, p*1e-2, tx11*1e-2, '/users/srio/Oasys/cubic_surface_step1.h5', overwrite=True)



#
#
#
#
# # step-2
# # to generate the working surface then use below shown line
#
# write_shadow_surface(Z1, p, -tx1, "/users/srio/Oasys/cubic_surface_step2.dat")
#
# # Divergence
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surface = ax.plot_surface(tx1, ty1, Z1, cmap='plasma')  # x,y,z
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# # ax.set_title('3D Surface Plot')
# # write_shadowSurface(Z1.T, Xc1,Yc1, "test_surf3.dat") # Z,x,y
# # Customize colorbar
# fig.colorbar(surface, ax=ax, shrink=0.5, aspect=10)
#
# # Show the plot
# plt.show()
# L = ty1[-1] - ty1[0]
#
# Lt = math.sqrt((ty1[-1] - ty1[0]) ** 2 + (tx1[-1] - tx1[0]) ** 2)
#
# D = 1200  # Distance between optic and source
#
# div2 = abs((L * np.sin(np.radians(45))) / D)  # These is divergence along the y axes of the optic
#
# div3 = (Lt * np.sin(np.radians(
#     45))) / D  # these is total divergence and given to the vertical divergence (-z) of the source and for vertical divergence (+z) is taken 0
#
# print("div2 = ", div2)
#
# print("div3 = ", div3)
#
# # print("div4 = ",np.degrees(div3 - div2))