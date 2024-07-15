from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import shutil
import os

def run_nvt_simulation(parm, inpcrd, pdb_file, temperature, nvt_steps, pdb_freq, cuda_device_index, nvt_output_pdb):

    # Convert temperature to kelvin
    temperature = temperature * kelvin
 
    # Initialize System and Integrator
    prmtop = AmberPrmtopFile(parm)
    inpcrd = AmberInpcrdFile(inpcrd)
    init_pdb = PDBFile(pdb_file)
    
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(temperature, 1/picosecond, 2*femtoseconds)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': cuda_device_index, 'CudaPrecision': 'mixed'}
    
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(init_pdb.positions)
    if inpcrd.boxVectors:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
    # Print initial energy
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy()
    print(f"Initial potential energy: {initial_energy}")
    
    # Minimize Energy
    simulation.minimizeEnergy(maxIterations=100000)
    
    # Print final energy
    state = simulation.context.getState(getEnergy=True)
    final_energy = state.getPotentialEnergy()
    print(f"Final potential energy: {final_energy}")
    
    # Set Reporters for NVT Simulation
    nvt_last_frame = nvt_output_pdb[:-4] + '_last_frame.pdb'
    simulation.reporters.append(PDBReporter(nvt_output_pdb, pdb_freq))
    simulation.reporters.append(PDBReporter(nvt_last_frame, nvt_steps))
    simulation.reporters.append(StateDataReporter(stdout, pdb_freq, step=True, time=True, potentialEnergy=True,
                                                  totalSteps=nvt_steps, temperature=True, progress=True, 
                                                  remainingTime=True, speed=True, separator='\t'))
    
    # Run NVT Simulation
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(nvt_steps)
    print("Finished NVT Simulation")

run_nvt_simulation(
    parm="parameterized_files/system.prmtop",
    inpcrd="parameterized_files/system.inpcrd",
    pdb_file="parameterized_files/system.pdb",
    temperature=300,
    nvt_steps=1000,
    pdb_freq=100,
    cuda_device_index='0',
    nvt_output_pdb="system_nvt_output.pdb")

# Create the equilibration directory if it does not exist
equilibration_dir = "simulation_openmm"
if not os.path.exists(equilibration_dir):
    os.makedirs(equilibration_dir)

# List of files to be moved
files_to_move = ["system_nvt_output.pdb", "system_nvt_output_last_frame.pdb"]

# Move the files to the equilibration directory
for file_name in files_to_move:
    if os.path.exists(file_name):
        shutil.move(file_name, os.path.join(equilibration_dir, file_name))
