from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Load the PDB file
fixer = PDBFixer(filename='/Users/edwardkim/IGEM/gh_scoresolver/tests/pdb/fabcomplexNOhydro.pdb')

# Find and add missing residues
fixer.findMissingResidues()

# Find and add missing atoms
fixer.findMissingAtoms()
fixer.addMissingHydrogens()

# Find non-standard residues and exclude UNL from replacement
fixer.findNonstandardResidues()
nonstandard_residues = fixer.nonstandardResidues
print("nonstards")
for residue_name in nonstandard_residues:
    print(residue_name)

if 'UNL' in nonstandard_residues:
    # Remove UNL from non-standard residues list
    del nonstandard_residues['UNL']


# Replace non-standard residues except UNL
fixer.replaceNonstandardResidues()

# Output the cleaned PDB file
with open('output_cleaned.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
