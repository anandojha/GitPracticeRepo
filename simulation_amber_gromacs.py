from openmm.unit import *
from openmm.app import *
from openmm import *
import parmed as pmd
import warnings
import shutil
import time
import sys
import os

def run_simulation(parm_file, inpcrd_file, pdb_file, temperature, nvt_steps, production_steps, cuda_device_index):

    # Function to read box vectors from the PDB file
    def read_box_vectors_from_pdb(pdb_filename):
        with open(pdb_filename, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("CRYST1"):
                    parts = line.split()
                    a, b, c = float(parts[1]), float(parts[2]), float(parts[3])
                    alpha, beta, gamma = float(parts[4]), float(parts[5]), float(parts[6])
                    return a, b, c, alpha, beta, gamma
        raise ValueError("CRYST1 record not found in PDB file")

    # Load Amber parameter and coordinate files
    parm = pmd.load_file(parm_file, inpcrd_file)
    # Save GROMACS topology and coordinate files
    parm.save('system.top', format='gromacs')
    parm.save('system.gro')

    # Read box vectors from the PDB file
    a, b, c, alpha, beta, gamma = read_box_vectors_from_pdb(pdb_file)
    # Define periodic box vectors based on CRYST1 record
    box_vectors = [Vec3(a, 0.0, 0.0), Vec3(0.0, b, 0.0), Vec3(0.0, 0.0, c)] * angstroms
    print(box_vectors)

    # Suppress specific warnings about consecutive residues with the same number
    warnings.filterwarnings("ignore", message="WARNING: two consecutive residues with same number")

    # Load the PDB file
    pdb = PDBFile(pdb_file)
    # Load the Gromacs topology file
    gromacs_topology = GromacsTopFile('system.top', periodicBoxVectors=box_vectors)

    # Create the OpenMM system
    system = gromacs_topology.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    # Integrator
    integrator = LangevinIntegrator(temperature*kelvin, 1/picosecond, 0.002*picoseconds)
    # Platform
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    # Create the Simulation object
    simulation = Simulation(gromacs_topology.topology, system, integrator, platform, properties)
    # Set the initial positions
    simulation.context.setPositions(pdb.positions)

    # Minimize the energy and print initial energy
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy()
    print(f'Initial potential energy: {initial_energy}')
    print('Minimizing...')
    simulation.minimizeEnergy(maxIterations=100000)
    # Print energy after minimization
    state = simulation.context.getState(getEnergy=True)
    minimized_energy = state.getPotentialEnergy()
    print(f'Potential energy after minimization: {minimized_energy}')

    # Equilibrate the system
    simulation.context.setVelocitiesToTemperature(temperature*kelvin)
    print('Equilibrating...')
    start_time = time.time()
    simulation.step(nvt_steps) 

    # Production run
    print('Running Production...')
    simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True, time=True, totalEnergy=True))
    simulation.reporters.append(DCDReporter('production.dcd', 1000))
    simulation.step(production_steps)
    end_time = time.time()
    print('Done!')

    # Print total simulation time
    total_time = end_time - start_time
    print(f'Total simulation time: {total_time} seconds')

    # Print final state
    state = simulation.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
    final_energy = state.getPotentialEnergy()
    kinetic_energy = state.getKineticEnergy()
    print(f'Final potential energy: {final_energy}')
    print(f'Final kinetic energy: {kinetic_energy}')

run_simulation(
    parm_file='parameterized_files/system.prmtop',
    inpcrd_file='parameterized_files/system.inpcrd',
    pdb_file='parameterized_files/system.pdb',
    temperature=300,
    nvt_steps=1000,
    production_steps=1000,
    cuda_device_index='0'
)

# Create the equilibration directory if it does not exist
equilibration_dir = "simulation_openmm_gromacs"
if not os.path.exists(equilibration_dir):
    os.makedirs(equilibration_dir)

# List of files to be moved
files_to_move = ["production.dcd", "system.gro", "system.top"]

# Move the files to the equilibration directory
for file_name in files_to_move:
    if os.path.exists(file_name):
        shutil.move(file_name, os.path.join(equilibration_dir, file_name))
