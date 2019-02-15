import numpy

wavelength = 0.7293
dspacing = 5.43/numpy.sqrt(3)

p = 39.0
R = 1.750

q = 1.0 / (1.0/p - 1.0/(-R))


thetaB = numpy.arcsin(wavelength/2/dspacing)
print("Bragg angle is: %f deg"%(thetaB*180/numpy.pi))



# sin2 theta1     sin2 theta2    sin theta1 + |sin theta2|
# ------------ + ------------- = -------------------------
#     p                q                  Rt

theta1 = (90-6.678633) * numpy.pi / 180
theta2 = (173.321367-90.0) * numpy.pi / 180

# theta1 = -6.678633 * numpy.pi / 180
# theta2 = 173.321367 * numpy.pi / 180

S1square = numpy.sin(theta1)**2
S2square = numpy.sin(theta2)**2
S1 = numpy.sin(theta1)
S2 = numpy.sin(theta2)




# q_shadow = 0.852
#
# f_inv =  S1square / p +  S2square / q_shadow
# Rt = (S1 + S2) / f_inv
#
# print("R nominal: ",R)
# print("q (lens): ",q)
# print("f(shadow): %f ; Rt (shadow): %f: "%(1.0/f_inv,Rt))

# q_shadow_calc =  S2square / ((S1+S2)/R - S2square/p )
q_shadow_calc =  S2square / ((S1+S2)/R + S2square/p )
print("q_shadow_calc",q_shadow_calc)



#
# krisch
#

# cos2 phi1     cos2 phi2    cos phi1 + |cos phi2|
# ------------ + ------------- = -------------------------
#     p                q                  Rt


phi1 = +6.678633 * numpy.pi / 180
phi2 = (180-173.321367) * numpy.pi / 180
C1square = numpy.cos(phi1)**2
C2square = numpy.cos(phi2)**2
C1 = numpy.cos(phi1)
C2 = numpy.cos(phi2)

q_krisch_calc =  C2square / ((C1+C2)/R + C2square/p )
print("q_krisch_calc",q_krisch_calc)

