import openmm as mm
from openmm import unit
import numpy as np
from sys import stdout
import matplotlib.pyplot as plt
import openmm.app as app
from openmm.app import Modeller
import subprocess
import os
import re
import sys
from progress.bar import IncrementalBar 
import pathlib

# & make new test folder
masterpath = os.path.join(pathlib.Path().absolute())
testpath = os.path.join(masterpath, "tests")
pdbname = "deca.pdb"
pdbpath = os.path.join(testpath,"pdb",pdbname)

if not os.path.exists(os.path.join(testpath, "pdb")):
    os.makedirs(os.path.join(testpath, "pdb"), exist_ok=True)

dir_list = os.listdir(testpath)

def increment_test_string(test_string):
    # regex to match folder names like 't{#}_{pdbname}'
    pattern = re.compile(r't(\d+)_' + re.escape(pdbname))
    
    match = pattern.match(test_string)
    if match:
        index = int(match.group(1))
        new_index = index + 1
        new_string = f"t{new_index}_{pdbname}"
        return new_string
    else:
        return None

def nextFolder(): 
    max_index = -1
    max_folder = None
    pattern = re.compile(r't(\d+)_' + re.escape(pdbname))

    for folder in dir_list:
        match = pattern.match(folder)
        if match:
            index = int(match.group(1))
            if index > max_index:
                max_index = index
                max_folder = folder

    if max_folder is None:
        return f"t0_{pdbname}"
    return increment_test_string(max_folder)

testname = nextFolder()
newpath = os.path.join(testpath, testname)

if not os.path.exists(newpath):
    os.makedirs(newpath, exist_ok=True)

os.makedirs(os.path.join(newpath, "cv"), exist_ok=True)
os.makedirs(os.path.join(newpath, "windows"), exist_ok=True)
os.makedirs(os.path.join(newpath, "hist"), exist_ok=True)

print(f"New folder created: '{newpath}'")

# & SYSTEM ####################################################################################################################################
pdb = app.PDBFile(pdbpath)

# ? if its PFOA or proteins of superset of PFOA, it might have some missing hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = app.ForceField('amber14-all.xml')
modeller.addHydrogens(forcefield)
pdb = modeller

# * system and simulation definition
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds, hydrogenMass=1.5*unit.amu)
integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.reporters.append(app.DCDReporter(os.path.join(newpath, 'hist', 'smd_traj.dcd'), 10000))
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.step(1000)

#&##################################################################################################################################

# ~ PARAMS ######################################################################################################

# * CV := distance betw CAs of two end residues
L_i = 1.3 # ? start length
L_f = 3.3 # ? end length
index1 = 8
index2 = 98
cv = mm.CustomBondForce('r')
cv.addBond(index1, index2)
num_win = 24

r0 = L_i*unit.nanometers #start value
fc_pull = 1000.0*unit.kilojoules_per_mole/unit.nanometers**2 #force constnat
v_pulling = 0.02*unit.nanometers/unit.picosecond #pulling velocity
dt = simulation.integrator.getStepSize() # simulation time step

total_steps = 30000 # ? total steps
increment_steps = 10 # ? step size
wTotal = 100000
wdelta = 1000 # ? record steps


# * Harmonic force definiton
pullingForce = mm.CustomCVForce('0.5 * fc_pull * (cv-r0)^2') #? maybe add rotational bias as well?
pullingForce.addGlobalParameter('fc_pull', fc_pull)
pullingForce.addGlobalParameter('r0', r0)
pullingForce.addCollectiveVariable("cv", cv)
system.addForce(pullingForce)
simulation.context.reinitialize(preserveState=True)

# * Window definition
windows = np.linspace(L_i, L_f, num_win)
window_coords = []
window_index = 0

print(f"Initialized System\nSystem: \n\t Forcefield Model: 'amber14-all.xml'\n \t Hydrogen Mass: '1.5 amu'\n \t Step size: {dt}\n\nParameters: \n\t {L_i} {L_f} {index1} {index2} {num_win}\n \t {fc_pull} {v_pulling}\n \t {total_steps} / {increment_steps} | {wTotal}/ {wdelta}")

n = index1
m = index2
assert n != m

with open(pdbpath, 'r') as pdb_file:
    atom_count = 0
    for line in pdb_file:
        if line.startswith("ATOM"):
            atom_count += 1
            if atom_count == n or atom_count == m:
                atom_number = int(line[6:11].strip())  # Atom serial number
                atom_name = line[12:16].strip()  # Atom name
                residue_name = line[17:20].strip()  # Residue name
                chain_id = line[21].strip()  # Chain identifier
                residue_number = int(line[22:26].strip())  # Residue sequence number
                x_coord = float(line[30:38].strip())  # X coordinate
                y_coord = float(line[38:46].strip())  # Y coordinate
                z_coord = float(line[46:54].strip())  # Z coordinate
                occupancy = float(line[54:60].strip())  # Occupancy
                temp_factor = float(line[60:66].strip())  # Temperature factor
                element_symbol = line[76:78].strip()  # Element symbol

                if atom_count == n :
                    print(f"index1: {index1}")
                else:
                    print(f"index2: {index2}")
                print(f"\t Atom Number: {atom_number}")
                print(f"\t Atom Name: {atom_name}")
                print(f"\t Residue Name: {residue_name}")
                print(f"\t Chain ID: {chain_id}")
                print(f"\t Residue Number: {residue_number}")
                print(f"\t Coordinates: ({x_coord}, {y_coord}, {z_coord})")
                print(f"\t Occupancy: {occupancy}")
                print(f"\t Temperature Factor: {temp_factor}")
                print(f"\t Element Symbol: {element_symbol}")
##~############################################################################################################################

#!########## SIMULATION ##################################################################################################################

# * SMD pulling loop
# Define the progress bar
progress_bar = IncrementalBar('Pull Loop', max=total_steps//increment_steps)

# Redirect stdout to a file or another stream
original_stdout = sys.stdout
with open('output.log', 'w') as f:
    sys.stdout = f

    # SMD pulling loop
    for i in range(total_steps//increment_steps):
        # Your simulation steps here
        simulation.step(increment_steps)
        current_cv_value = pullingForce.getCollectiveVariableValues(simulation.context)
        r0 += v_pulling * dt * increment_steps
        simulation.context.setParameter('r0', r0)
        if (window_index < len(windows) and current_cv_value >= windows[window_index]):
            window_coords.append(simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions())
            window_index += 1
        progress_bar.next()

    # Finish the progress bar
    progress_bar.finish()

# Reset stdout to its original value
sys.stdout = original_stdout

#simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, time=True, potentialEnergy=True, temperature=True, speed=True))

# * Save Windows
for i, coords in enumerate(window_coords):
   # Save window structures to 'windows' directory
   window_outfile = open(os.path.join(newpath,"windows", f'window_{i}.pdb'), 'w')
   app.PDBFile.writeFile(simulation.topology, coords, window_outfile)
   window_outfile.close()

#!##################################################################################################################################

##^## RUN WINDOWS AND REUSE SIMULATION -> CV TIME SERIES FILES ########################################################################################################################################

def run_window(window_index):
    pdb = app.PDBFile(os.path.join(newpath, "windows", f'window_{window_index}.pdb'))  # Load from 'windows' directory
    simulation.context.setPositions(pdb.positions)
    r0 = windows[window_index]
    simulation.context.setParameter('r0', r0)
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    simulation.step(wdelta)
    total_steps = wTotal
    record_steps = wdelta
    cv_values = []

    # Define the progress bar
    progress_bar = IncrementalBar(f'Window {window_index}', max=total_steps//record_steps)

    # SMD pulling loop
    for i in range(total_steps//record_steps):
        simulation.step(record_steps)
        current_cv_value = pullingForce.getCollectiveVariableValues(simulation.context)
        cv_values.append([i, current_cv_value[0]])
        
        # Update the progress bar
        progress_bar.next()

    # Finish the progress bar
    progress_bar.finish()

    # Save CV values to file
    np.savetxt(os.path.join(newpath, "cv", f'cv_{window_index}.txt'), np.array(cv_values))  # Save to 'cv' directory


for n in range(num_win):
   run_window(n)

#^#######################################################################################################################################

# & histograms write and plot ##########################################################################################################################################
metafilelines = []
for i in range(len(windows)):
    cvpath = os.path.join(newpath,"cv",f'cv_{i}.txt')
    data = np.loadtxt(f'{cvpath}')
    plt.hist(data[:,1])
    metafileline = f'{cvpath} {windows[i]} {wdelta} \n' # ? do i do just 'cv_{i}.txt {windows[i]}' or '{cvpath}' ?
    metafilelines.append(metafileline)

plt.title("Histogram")
plt.xlabel("r (nm)")
plt.ylabel("count")

with open(os.path.join(newpath,"hist","metafile.txt"), "w") as f:
    f.writelines(metafilelines)
output_directory = os.pathjoin(newpath,"results")
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
plt.savefig(output_directory)
#&#######################################################################################################################################################################################

# ~ WHAM on metafile.txt ################################################################################################################################################################

# ! USE POPEN TO EXECUTE child program ---->>>>> /Users/edwardkim/Downloads/wham 1.3 3.3 50 1e-6 300 0 metafile.txt pmf.txt > wham_log.txt

wham_xecpath = os.path.join(masterpath,"wham4","wham","wham")
print([wham_xecpath, str(L_i), str(L_f), '50', '1e-6', '300', '0', os.path.join(newpath,"hist", "metafile.txt"), os.path.join(newpath,"hist", "pmf.txt")])
args = [wham_xecpath, str(L_i), str(L_f), '50', '1e-6', '300', '0', os.path.join(newpath,"hist", "metafile.txt"), os.path.join(newpath,"hist", "pmf.txt")]

# give perms to the path folder 
os.chmod(wham_xecpath, 0o755)
# process obj
process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# wait for process termination and capture output
stdout_data, stderr_data = process.communicate()

# write output to wham logs
with open(os.path.join(newpath,"hist","wham_log.txt"), "w") as output_file:
    output_file.write(stdout_data.decode())  # Decode stdout_data from bytes to string

# check for errors
if stderr_data:
    print("Error occurred:", stderr_data.decode())  # Decode stderr_data from bytes to string

# ~ ###############################################################################################################################################################

# & plot PMF

coordinates = []
free_energies = []
probabilities = []

with open(os.path.join(newpath,"hist", "pmf.txt"), 'r') as file:
    for line in file:
        # skip comment lines ("# ...")
        if line.startswith('#'):
            continue
        # split lines into columns
        columns = line.split()
        if len(columns) == 5:
            coor = float(columns[0])
            free = float(columns[1]) if columns[1] != 'nan' else None
            prob = float(columns[3]) if columns[3] != 'nan' else None
            coordinates.append(coor)
            free_energies.append(free)
            probabilities.append(prob)

fig, ax1 = plt.subplots()

ax1.set_xlabel('ξ') 
ax1.set_ylabel('Free Energy', color='tab:blue')
ax1.plot(coordinates, free_energies, 'b-', label='Free Energy')
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax2 = ax1.twinx()
ax2.set_ylabel('Probability', color='tab:red')
ax2.plot(coordinates, probabilities, 'r-', label='Probability')
ax2.tick_params(axis='y', labelcolor='tab:red')

fig.tight_layout()
fig.legend(loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)

plt.title('Histogram Data Plot')
plt.show()


plt.savefig(output_directory)
