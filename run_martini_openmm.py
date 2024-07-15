import os

def run_command(command):
    """
    Run a system command and print it.

    Parameters:
    command (str): The command to execute.
    """
    print(f"Executing: {command}")
    os.system(command)

# Copy the files to the martini simulation directory
current_dir = os.getcwd()

run_command("cp -r " + current_dir + "/" + "martini_openmm.py" + " " +  current_dir + "/" + "openmm_martini_sims")

# Run OpenMM simulations
os.chdir(current_dir + "/" + "openmm_martini_sims")

# Check if the martini_openmm directory already exists, delete if it does
if os.path.exists("martini_openmm"):
    run_command("rm -rf martini_openmm")

# Clone the repository
run_command("git clone https://github.com/maccallumlab/martini_openmm.git")

# Run the simulation script
run_command("python martini_openmm.py")

# Return to the original directory
os.chdir(current_dir)


