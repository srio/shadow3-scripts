import numpy
from srxraylib.plot.gol import plot, plot_image
from srxraylib.metrology.dabam import write_shadowSurface


filein = "/users/srio/Oasys/diaboloid_correction.txt"

sagittal = numpy.loadtxt(filein)
sagittal_x = sagittal[:,0].copy()
sagittal_y = sagittal[:,1].copy() * 1

# print(sagittal_x.shape)
plot(sagittal_x, sagittal_y)

tangential_x = numpy.linspace(-0.4, 0.4, 1001)
tangential_y = numpy.zeros_like(tangential_x)
mesh = numpy.zeros((sagittal_x.size, tangential_x.size))

for i in range(tangential_x.size):
    mesh[:,i] = sagittal_y

mesh -= mesh[mesh.shape[0]//2, mesh.shape[1]//2]
plot_image(mesh, sagittal_x, tangential_x, aspect='auto')


fileout = "/users/srio/Oasys/diaboloid_bl1222_goldenrule_shadow.dat"
print(mesh.shape, sagittal_x.shape, tangential_x.shape)
write_shadowSurface(mesh.T, sagittal_x, tangential_x, outFile=fileout)

mirror_size = sagittal_x[-1] * 2
d_mirror_over_crystal = 18.8 / 14.1
crystal_size =  mirror_size / d_mirror_over_crystal
print("mirror_size, crystal size: ", mirror_size, crystal_size)


fileout = "/users/srio/Oasys/diaboloid_crystal.dat"
print(mesh.shape, sagittal_x.shape, tangential_x.shape)
bragg_angle = (90 - 80.824081) * numpy.pi / 180
correction_factor = (0.002 / bragg_angle) * 1.0
print("Correction factor: ", correction_factor)
write_shadowSurface(mesh.T * correction_factor, sagittal_x / d_mirror_over_crystal, tangential_x / 0.4 * 0.015, outFile=fileout)