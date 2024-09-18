# Script to get atom number of selected atom in ChimeraX

# Import necessary modules
from chimerax.core.commands import run
from chimerax.atomic import selected_atoms


# Get the selected atom
atoms = selected_atoms(session)
if atoms:
    selected_atom = atoms[0]  # Assuming only one atom is selected
    # Print the atom number (serial number)
    print(f"Atom number (serial number): {selected_atom.serial_number}")
else:
    print("No atom selected.")

# runscript /Users/edwardkim/IGEM/gh_scoresolver/atomnumber_chimerax.py