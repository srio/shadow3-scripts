import numpy
import xraylib
energy = 15.0


delta_Be = 1.0 - xraylib.Refractive_Index_Re("Be", energy, 1.848)
delta_diamond = 1.0 - xraylib.Refractive_Index_Re("C", energy, 3.51)
delta_Al = 1.0 - xraylib.Refractive_Index_Re("Al", energy, 2.7)

print("delta Be", delta_Be)
print("delta Diamond", delta_diamond)
print("delta Al", delta_Al)

f = 0.627
radii = [200e-6, 100e-6, 50e-6, 30e-6]

for radius in radii:
    print("================ radius[um] = ", 1e6 * radius)
    print("Be lenses: ", int(radius / (2 * delta_Be * f)))
    print("Diamond lenses: ", int(radius / (2 * delta_diamond * f)))
    print("Al lenses: ", int(radius / (2 * delta_Al * f)))
