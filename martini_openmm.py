import time
import sys
import os

# Ensure the inner martini_openmm directory is in the Python path
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, "martini_openmm"))

# OpenMM simulations in Martini
# URL: https://github.com/maccallumlab/martini_openmm/tree/master

# Import necessary modules
from mdtraj.reporters import XTCReporter
import martini_openmm as martini
from openmm.unit import *
from openmm.app import *
from sys import stdout
from openmm import *

def run_openMM(epsilon_r):
    """
    Run the OpenMM simulation with Martini force field.

    Parameters:
    epsilon_r (float): Relative permittivity for the simulation.
    """
    platform = Platform.getPlatformByName("CUDA")
    properties = {'Precision': 'double'}
    conf = GromacsGroFile("solvated.gro")
    box_vectors = conf.getPeriodicBoxVectors()
    defines = {}
    try:
        with open("defines.txt") as def_file:
            for line in def_file:
                line = line.strip()
                defines[line] = True
    except FileNotFoundError:
        pass

    top = martini.MartiniTopFile("system.top", periodicBoxVectors=box_vectors, defines=defines, epsilon_r=epsilon_r)
    system = top.create_system(nonbonded_cutoff=1.1 * nanometer)
    integrator = LangevinIntegrator(300 * kelvin, 10.0 / picosecond, 30 * femtosecond)
    integrator.setRandomNumberSeed(0)
    simulation = Simulation(top.topology, system, integrator, platform, properties)
    simulation.context.setPositions(conf.getPositions())
    ### Minimization ###
    simulation.reporters.append(PDBReporter('minimization.pdb', 100))
    simulation.reporters.append(StateDataReporter(stdout, 500, step=True, potentialEnergy=True, temperature=True, volume=True))
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=50000, tolerance=1.0)
    energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print("System minimized at", energies, "\n")
    ### NVT Equilibration ###
    simulation.context.setVelocitiesToTemperature(300 * kelvin)
    print('Running NVT equilibration...')
    simulation.step(50000)
    ### NPT Equilibration ###
    system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin))
    simulation.context.reinitialize(True)
    print('Running NPT equilibration...')
    simulation.step(50000)
    simulation.saveState('equilibration.state')
    simulation.saveCheckpoint('equilibration.chk')
    ### Production Run ###
    simulation.reporters.append(StateDataReporter("production.log", 1000, step=True, potentialEnergy=True, totalEnergy=True, density=True, temperature=True, volume=True))
    xtc_reporter = XTCReporter('production.xtc', 500)
    simulation.reporters.append(xtc_reporter)
    # Measure the performance of the production run
    production_steps = 100000  # Number of steps for the production run
    timestep = 30 * femtoseconds  # Timestep size
    simulated_time = production_steps * timestep  # Total simulated time
    print("Running simulation...")
    start_time = time.time()  # Start timing
    simulation.step(production_steps)
    end_time = time.time()  # End timing
    # Calculate the performance
    elapsed_time = end_time - start_time  # Elapsed wall-clock time in seconds
    simulated_ns = simulated_time.value_in_unit(nanoseconds)  # Simulated time in nanoseconds
    performance_ns_per_day = simulated_ns / (elapsed_time / 86400)  # Performance in ns/day
    print(f"Production run took {elapsed_time:.2f} seconds")
    print(f"Simulated time: {simulated_ns:.2f} ns")
    print(f"Performance: {performance_ns_per_day:.2f} ns/day")

# Execute the run function with a relative permittivity of 15
run_openMM(15)

