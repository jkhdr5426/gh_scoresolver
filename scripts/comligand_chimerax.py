from chimerax.core.commands import run
from chimerax.atomic import AtomicStructure
import numpy as np

def calculate_center_of_mass(atoms):
    coordinates = np.array([atom.scene_coord for atom in atoms])
    masses = np.array([atom.element.mass for atom in atoms])
    center_of_mass = np.average(coordinates, axis=0, weights=masses)
    return center_of_mass

# Get the current structure
structures = list(session.models.list(type=AtomicStructure))
if not structures:
    raise ValueError("No structure loaded")
structure = structures[0]

# Select the ligand
# run(session, "select ::name='UNK'")

# Get the selected atoms (the ligand atoms)
selected_atoms = structure.atoms[structure.atoms.selected]

# Calculate the center of mass of the ligand
if len(selected_atoms) == 0:
    raise ValueError("No ligand selected")
center_of_mass = calculate_center_of_mass(selected_atoms)

# Find the atom closest to the center of mass
closest_atom = None
min_distance = float('inf')
for atom in structure.atoms:
    distance = np.linalg.norm(atom.scene_coord - center_of_mass)
    if distance < min_distance:
        min_distance = distance
        closest_atom = atom

# Print the result
if closest_atom:
    atom_number = closest_atom.serial_number
    atom_name = closest_atom.name
    atom_element = closest_atom.element.name
    residue = closest_atom.residue
    print(f"Closest atom: {atom_name} ({atom_element}) in residue {residue.name} {residue.number}, atom number: {atom_number}, at distance {min_distance}")
    run(session, f"select {closest_atom.atomspec}")
else:
    print("No atoms found in the structure")

# ! type in chimera:
# runscript <copy absolute path of this script>