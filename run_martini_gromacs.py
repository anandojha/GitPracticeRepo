from simtk.openmm.app import PDBFile
from pdbfixer import PDBFixer
import subprocess
import shutil
import os

def run_command(command):
    """
    Run a system command and print it.

    Parameters:
    command (str): The command to execute.
    """
    print(f"Executing: {command}")
    os.system(command)

def count_water_molecules(gro_file):
    """
    Count the number of water molecules in a GRO file.

    Parameters:
    gro_file (str): Path to the GRO file.

    Returns:
    int: The number of water molecules.
    """
    water_count = 0
    with open(gro_file, 'r') as file:
        for line in file:
            if " W" in line:
                water_count += 1
    return water_count

def remove_HOH(input_file, output_file):
    """
    Remove water molecules (HOH) from a PDB file.

    Parameters:
    input_file (str): Path to the input PDB file.
    output_file (str): Path to the output PDB file with HOH removed.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as file:
        for line in lines:
            if "HOH" not in line:
                file.write(line)

def remove_NI(input_file, output_file):
    """
    Remove nickel ions (NI) from a PDB file.

    Parameters:
    input_file (str): Path to the input PDB file.
    output_file (str): Path to the output PDB file with NI removed.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as file:
        for line in lines:
            if " NI " not in line:
                file.write(line)

def remove_CRYST(input_file, output_file):
    """
    Remove CRYST1 lines from a PDB file.

    Parameters:
    input_file (str): Path to the input PDB file.
    output_file (str): Path to the output PDB file with CRYST1 removed.
    """
    with open(input_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as file:
        for line in lines:
            if "CRYST1" not in line:
                file.write(line)

def update_topology_file(top_file_path, num_water_molecules):
    """
    Update the topology file with the number of water molecules.

    Parameters:
    top_file_path (str): Path to the topology file.
    num_water_molecules (int): The number of water molecules to add to the topology.
    
    Returns:
    str: Path to the updated topology file.
    """
    # Read the topology file
    with open(top_file_path, 'r') as file:
        top_content = file.readlines()
    
    # Find the [ molecules ] section
    molecules_section_found = False
    updated_content = []
    for line in top_content:
        updated_content.append(line)
        if "[ molecules ]" in line:
            molecules_section_found = True

    # Append the water molecule line correctly after the last entry
    if molecules_section_found:
        for i in range(len(updated_content) - 1, -1, -1):
            if updated_content[i].strip() and not updated_content[i].startswith(";") and not updated_content[i].startswith("W"):
                updated_content.insert(i + 1, f"\nW             {num_water_molecules}\n")
                break

    # Write the updated content back to the original file to overwrite it
    with open(top_file_path, 'w') as file:
        file.writelines(updated_content)
    
    return top_file_path

def process_pdb(input_pdb_path, output_pdb_path):
    """
    Process a PDB file to remove existing TER lines and add TER lines between chains A and B, and B and I.

    Parameters:
    input_pdb_path (str): Path to the input PDB file.
    output_pdb_path (str): Path to the output modified PDB file.
    """
    # Read the PDB file
    with open(input_pdb_path, 'r') as file:
        pdb_lines = file.readlines()
    
    # Process the PDB file to remove existing TER lines and add TER lines between chains A and B, and B and I
    modified_pdb_lines = []
    current_chain = None
    for line in pdb_lines:
        if line.startswith("TER"):
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain_id = line[21]
            if current_chain is None:
                current_chain = chain_id
            elif current_chain != chain_id:
                if (current_chain == 'A' and chain_id == 'B') or (current_chain == 'B' and chain_id == 'I'):
                    modified_pdb_lines.append("TER\n")
                current_chain = chain_id
        modified_pdb_lines.append(line)
    
    # Ensure the final TER is added before the END line
    if modified_pdb_lines[-1].startswith("END"):
        modified_pdb_lines.insert(-1, "TER\n")
    else:
        modified_pdb_lines.append("TER\n")
    
    # Write the modified PDB file
    with open(output_pdb_path, 'w') as file:
        file.writelines(modified_pdb_lines)

def fix_missing_atoms(input_pdb, output_pdb):
    """
    Check and fix missing atoms in a PDB file using PDBFixer.

    Parameters:
    input_pdb (str): Path to the input PDB file.
    output_pdb (str): Path to the output PDB file with fixed atoms.
    """
    # Initialize PDBFixer with the input PDB file
    fixer = PDBFixer(filename=input_pdb)

    # Find and add missing residues
    fixer.findMissingResidues()

    # Find and add missing atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Write the fixed PDB to the output file
    with open(output_pdb, 'w') as output_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, output_file)

    print(f"Missing atoms fixed and output written to {output_pdb}")
    
# Define the directories to delete (if existing) and recreate
directories = ['production', 'minimization', 'equilibration', 'intermediate', 'openmm_martini_sims']

# Delete the specified directories if they exist
for directory in directories:
    if os.path.exists(directory):
        shutil.rmtree(directory)
        print(f"Deleted {directory} directory")

# Recreate the specified directories
for directory in directories:
    os.makedirs(directory, exist_ok=True)
    print(f"Created {directory} directory")
    
# Copy the input files to the current directory
current_dir = os.getcwd()
run_command("cp -r " + current_dir + "/" + "martini_input_files" + "/* " +  current_dir)

# Step 0: Download the PDB file, decompress it, and rename it
run_command("wget https://files.rcsb.org/download/2XWJ.pdb1.gz")
run_command("gunzip 2XWJ.pdb1.gz")
run_command("mv 2XWJ.pdb1 2XWJ.pdb")

# Step 1: Copy the PDB file
# This command copies the original PDB file to a new file for processing
run_command("cp -r 2XWJ.pdb system_I.pdb")

# Step 2: Process the PDB file with pdb4amber
# This command processes the PDB file using pdb4amber, removes unnecessary connections and reduces the file
run_command("pdb4amber -i system_I.pdb -o system_II.pdb --no-conect --reduce -s :NAG")
# Remove temporary and unnecessary files generated during pdb4amber processing
run_command("rm -rf *nonprot* *renum* *sslink* *log*")

# Step 3: Remove water molecules, NI ions and crystal information from the PDB file using Python functions
# This script removes water molecules from the PDB file
remove_HOH('system_II.pdb', 'system_III.pdb')
remove_NI('system_III.pdb', 'system_IV.pdb')
remove_CRYST('system_IV.pdb', 'system_V.pdb')

# Step 4: Process the PDB file to remove existing TER lines and add TER lines between chains
process_pdb(input_pdb_path="system_V.pdb", output_pdb_path="system_VI.pdb")

# Step 5: Check and fix missing atoms using PDBFixer
fix_missing_atoms(input_pdb="system_VI.pdb", output_pdb="system.pdb")

# Run martinize.py to Generate Structure and Topology Files
# This command uses martinize.py to convert the PDB file into a coarse-grained model using the Martini force field
run_command("python martinize.py -f system.pdb -o system.top -x system-CG.pdb -p backbone -ff martini22")

# Generate a box with gmx editconf
# This command centers the coarse-grained model in a cubic box with a minimum distance of 0.5 nm to the box edges
run_command("gmx editconf -f system-CG.pdb -o system-CG.gro -c -d 0.10 -bt cubic")

# Prepare the input files for minimization
# This command prepares the GROMACS input file for energy minimization
run_command("gmx grompp -p system.top -f minimization.mdp -c system-CG.gro -o minimization-vac.tpr")

# Run the minimization
# This command runs the energy minimization in GROMACS
run_command("gmx mdrun -deffnm minimization-vac -v")

# Ensure that the water box file is available and solvate the system with gmx solvate
# This command solvates the minimized system by adding water molecules
run_command("gmx solvate -cp minimization-vac.gro -cs water.gro -o solvated.gro")

# Count the number of water molecules in the solvated system
num_water_molecules = count_water_molecules("solvated.gro")
print(num_water_molecules)

# Update the topology file with the number of water molecules
# This updates the topology file to include the correct number of water molecules
updated_top_file_path = update_topology_file(top_file_path="system.top", num_water_molecules=num_water_molecules)
print(f"Updated topology file saved to: {updated_top_file_path}")

# Prepare the input files for minimization again, now with the solvated system
run_command("gmx grompp -p system.top -c solvated.gro -f minimization.mdp -o minimization.tpr")

# Run the minimization on the solvated system
run_command("gmx mdrun -deffnm minimization -v")

# Prepare the input files for equilibration with the given mdp file
# This command prepares the GROMACS input file for equilibration
run_command("gmx grompp -p system.top -c minimization.gro -f equilibration.mdp -o equilibration.tpr -maxwarn 2")

# Run the equilibration
# This command runs the equilibration in GROMACS
run_command("gmx mdrun -deffnm equilibration -v")

# Prepare the input files for the production run
# This command prepares the GROMACS input file for the production run
run_command("gmx grompp -p system.top -c equilibration.gro -f dynamic.mdp -o dynamic.tpr -maxwarn 1")

# Run the production simulation
# This command runs the production molecular dynamics simulation in GROMACS
run_command("gmx mdrun -deffnm dynamic -v")

# Create an index file excluding water
make_ndx_command = "echo 'del 0\nkeep 0\nq\n' | gmx make_ndx -f dynamic.gro -o dynamic.ndx"

# Convert the trajectory to PDB focusing on the solute group
trjconv_command = "echo '0' | gmx trjconv -s dynamic.tpr -f dynamic.xtc -o dynamic.pdb -n dynamic.ndx -pbc mol -center"

# Define the commands to run
commands = [make_ndx_command,trjconv_command]

# Run the commands
for cmd in commands:
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(process.stdout)
    print(process.stderr)

# Delete unnecessary files
run_command("rm -rf *#* *prev*")

# Create directories for minimization, equilibration, production and intermediate files
os.makedirs('production', exist_ok=True)
os.makedirs('minimization', exist_ok=True)
os.makedirs('equilibration', exist_ok=True)
os.makedirs('intermediate', exist_ok=True)
os.makedirs('openmm_martini_sims', exist_ok=True)

# Define the files to move for each category
production_files = [
    'dynamic.cpt', 'dynamic.edr', 'dynamic.gro', 'dynamic.log', 
    'dynamic.mdp', 'dynamic.ndx', 'dynamic.pdb', 'dynamic.tpr', 
    'dynamic.xtc']

minimization_files = [
    'minimization.edr', 'minimization.gro', 'minimization.log', 
    'minimization.mdp', 'minimization.tpr', 'minimization.trr', 
    'minimization-vac.edr', 'minimization-vac.gro', 'minimization-vac.log', 
    'minimization-vac.tpr', 'minimization-vac.trr']

equilibration_files = [
    'equilibration.cpt', 'equilibration.edr', 'equilibration.gro', 
    'equilibration.log', 'equilibration.mdp', 'equilibration.tpr', 
    'equilibration.xtc']

intermediate_files = [
    'system_I.pdb', 'system_II.pdb', 'system_III.pdb', 
    'system_IV.pdb', 'system_V.pdb', 'system_VI.pdb',
    'system.pdb', '2XWJ.pdb', 'system-CG.pdb', 'system-CG.gro']

openmm_martini_files = ['Protein_A.itp','Protein_B.itp',"martini.itp",
                        'Protein_C.itp','solvated.gro','system.top']

# Move production files
for file in production_files:
    if os.path.exists(file):
        shutil.move(file, 'production')
    else:
        print(f"{file} does not exist")

# Move minimization files
for file in minimization_files:
    if os.path.exists(file):
        shutil.move(file, 'minimization')
    else:
        print(f"{file} does not exist")

# Move equilibration files
for file in equilibration_files:
    if os.path.exists(file):
        shutil.move(file, 'equilibration')
    else:
        print(f"{file} does not exist")
        
# Move intermediate_files
for file in intermediate_files:
    if os.path.exists(file):
        shutil.move(file, 'intermediate')
        print(f"Moved {file} to intermediate")
    else:
        print(f"{file} does not exist")
        
# Move openmm_martini files
for file in openmm_martini_files:
    if os.path.exists(file):
        shutil.move(file, 'openmm_martini_sims')
        print(f"Moved {file} to openmm_martini_sims")
    else:
        print(f"{file} does not exist")

print("Files have been moved to their respective directories.")

# List of files to delete
files_to_delete = ['martinize.py', 'water.gro', 'mdout.mdp']

# Delete the specified files
for file in files_to_delete:
    if os.path.exists(file):
        os.remove(file)
        print(f"Deleted {file}")
    else:
        print(f"{file} does not exist")

print("Specified files have been deleted.")

