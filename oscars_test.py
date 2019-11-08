
# coding: utf-8

# Plots inline for notebook
#get_ipython().run_line_magic('matplotlib', 'inline')

# Import the OSCARS SR module
from srxraylib.plot.gol import set_qt
import oscars.sr
from oscars.plots_mpl import *
from oscars.parametric_surfaces import PSCylinder

import numpy

def undulator_spectrum():
    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear all existing fields and create an undulator field
    osr.clear_bfields()
    osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)

    # Define simple electron beam
    osr.set_particle_beam(energy_GeV=3, x0=[0, 0, -1], current=0.5)

    # Define the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # Calculate spectrum at 30 [m]
    spectrum = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[100, 2000])

    # Optionally import the plotting tools (matplotlib)


    # Plot spectrum
    plot_spectrum(spectrum)

def undulator_flux():
    # # coding: utf-8
    #
    # # Plots inline for notebook
    # # get_ipython().run_line_magic('matplotlib', 'inline')
    #
    # # Import the OSCARS SR module
    # import oscars.sr
    #
    # # Optionally import the plotting tools (matplotlib)
    # from oscars.plots_mpl import *

    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear all existing fields and create an undulator field
    osr.clear_bfields()
    osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)

    # Define simple electron beam
    osr.set_particle_beam(energy_GeV=3, x0=[0, 0, -1], current=0.5)

    # Define the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # Calculate spectrum at 30 [m].  Note use of the nthreads argument.
    flux = osr.calculate_flux_rectangle(
        plane='XY',
        energy_eV=143.8,
        width=[0.01, 0.01],
        npoints=[101, 101],
        translation=[0, 0, 30]
    )

    # Plot flux
    plot_flux(flux)

def undulator_power_density():
    # # coding: utf-8
    #
    # # Plots inline for notebook
    # get_ipython().run_line_magic('matplotlib', 'inline')
    #
    # # Import the OSCARS SR module
    # import oscars.sr
    #
    # # Optionally import the plotting tools (matplotlib)
    # from oscars.plots_mpl import *

    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear all existing fields and create an undulator field
    osr.clear_bfields()
    osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)

    # Define simple electron beam
    osr.set_particle_beam(energy_GeV=3, x0=[0, 0, -1], current=0.5)

    # Define the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # Calculate spectrum at 30 [m].  Note use of the nthreads argument.
    power_density = osr.calculate_power_density_rectangle(
        plane='XY',
        width=[0.05, 0.05],
        npoints=[101, 101],
        translation=[0, 0, 30]
    )

    # Plot power density
    plot_power_density(power_density)

def undulator_3d_power_density():
    # # coding: utf-8
    #
    # # Plots inline for notebook
    # get_ipython().run_line_magic('matplotlib', 'inline')
    #
    # # Import the OSCARS SR module
    # import oscars.sr
    #
    # # Import the 3D and parametric surfaces utilities
    # from oscars.plots3d_mpl import *
    # from oscars.parametric_surfaces import *

    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear all existing fields and create an undulator field
    osr.clear_bfields()
    osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.050], nperiods=31)

    # Define simple electron beam
    osr.set_particle_beam(energy_GeV=3, x0=[0, 0, -1], current=0.5)

    # Define the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # First create the surface of interest
    cylinder = PSCylinder(R=0.020, L=0.010, nu=101, nv=101)

    # Run calculation and plotting
    pd = power_density_3d(osr, cylinder, rotations=[osr.pi() / 2, 0, 0], translation=[0, 0, 30])

def example_032_undulator_flux():
    # # Import the OSCARS SR module
    # import oscars.sr
    #
    # # Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
    # from oscars.plots_mpl import *

    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear any existing fields (just good habit in notebook style) and add an undulator field
    osr.clear_bfields()
    osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

    # Just to check the field that we added seems visually correct
    plot_bfield(osr)

    # Setup beam similar to NSLSII
    osr.clear_particle_beams()
    osr.set_particle_beam(x0=[0, 0, -1], energy_GeV=3, current=0.500)

    # Set the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # Run the particle trajectory calculation
    trajectory = osr.calculate_trajectory()

    # Plot the trajectory position and velocity
    plot_trajectory_position(trajectory)
    plot_trajectory_velocity(trajectory)

    # Calculate spectrum zoom
    spectrum = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[145, 160], npoints=200)
    plot_spectrum(spectrum)

    flux = osr.calculate_flux_rectangle(
        plane='XY',
        energy_eV=153,
        width=[0.01, 0.01],
        npoints=[101, 101],
        translation=[0, 0, 30]
    )

    plot_flux(flux)


def example_042_undulator_power_density():
    # # Import the OSCARS SR module
    # import oscars.sr
    #
    # # Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
    # from oscars.plots_mpl import *

    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear any existing fields (just good habit in notebook style) and add an undulator field
    osr.clear_bfields()
    osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

    # Just to check the field that we added seems visually correct
    plot_bfield(osr)

    # Setup beam similar to NSLSII
    osr.clear_particle_beams()
    osr.set_particle_beam(x0=[0, 0, -1], energy_GeV=3, current=0.500)

    # Set the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # Run the particle trajectory calculation
    trajectory = osr.calculate_trajectory()

    # Plot the trajectory position and velocity
    plot_trajectory_position(trajectory)
    plot_trajectory_velocity(trajectory)

    power_density = osr.calculate_power_density_rectangle(
        plane='XY',
        width=[0.05, 0.05],
        npoints=[101, 101],
        translation=[0, 0, 30]
    )

    plot_power_density(power_density)

def example_001_dipole_trajectory():
    # Import the OSCARS SR module
    # import oscars.sr
    #
    # # Import basic plot utilities.  You don't need these to run OSCARS, but it's used here for basic plots
    # from oscars.plots_mpl import *

    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear any existing fields (just good habit in notebook style) and add an undulator field
    osr.clear_bfields()
    osr.add_bfield_uniform(bfield=[0, -0.4, 0], width=[0, 0, 1])

    # Just to check the field that we added seems visually correct
    plot_bfield(osr)

    # Setup beam similar to NSLSII
    osr.clear_particle_beams()
    osr.set_particle_beam(x0=[0, 0, -1], energy_GeV=3, current=0.500)

    # Set the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # Verify input information - print all to screen
    osr.print_all()

    # Run the particle trajectory calculation
    trajectory = osr.calculate_trajectory()

    # Plot the trajectory position and velocity
    plot_trajectory_position(trajectory)
    plot_trajectory_velocity(trajectory)

    # Setup beam similar to NSLSII
    osr.clear_particle_beams()
    osr.set_particle_beam(energy_GeV=3, current=0.500)

    # Set the start and stop times for the calculation
    osr.set_ctstartstop(-1, 1)

    # Run the particle trajectory calculation
    trajectory = osr.calculate_trajectory()

    # Plot the trajectory position and velocity
    plot_trajectory_position(trajectory)
    plot_trajectory_velocity(trajectory)

def undulator_radiation_srio():
    # # Import the OSCARS SR module
    # import oscars.sr
    #
    # # Import basic plot utilities (matplotlib).  You don't need these to run OSCARS, but it's used here for basic plots
    # from oscars.plots_mpl import *

    # Create a new OSCARS object.  Default to 8 threads and always use the GPU if available
    osr = oscars.sr.sr(nthreads=8, gpu=1)

    # Clear any existing fields (just good habit in notebook style) and add an undulator field
    osr.clear_bfields()
    osr.add_bfield_undulator(bfield=[0, 1, 0], period=[0, 0, 0.049], nperiods=21)

    # Just to check the field that we added seems visually correct
    plot_bfield(osr)

    # Setup beam similar to NSLSII
    osr.clear_particle_beams()
    osr.set_particle_beam(x0=[0, 0, -1], energy_GeV=3, current=0.500)

    # Set the start and stop times for the calculation
    osr.set_ctstartstop(0, 2)

    # Run the particle trajectory calculation
    trajectory = osr.calculate_trajectory()

    # Plot the trajectory position and velocity
    plot_trajectory_position(trajectory)
    plot_trajectory_velocity(trajectory)

    # Calculate spectrum zoom
    spectrum = osr.calculate_spectrum(obs=[0, 0, 30], energy_range_eV=[145, 160], npoints=200)
    # print(">>>",spectrum)
    plot_spectrum(spectrum)

    flux = osr.calculate_flux_rectangle(
        plane='XY',
        energy_eV=153,
        width=[0.01, 0.01],
        npoints=[101, 101],
        translation=[0, 0, 30]
    )

    plot_flux(flux)
    print(">>>", flux,type(flux))

    print(numpy.array(flux).shape)
if __name__ == "__main__":

    set_qt()

    # undulator_spectrum()
    # undulator_flux()
    # undulator_power_density()
    # undulator_3d_power_density()   #??????????
    # example_032_undulator_flux()
    # example_042_undulator_power_density()
    # example_001_dipole_trajectory()

    undulator_radiation_srio()