from openmm import app
from openmm import XmlSerializer

# Load the ligand parameters into OpenMM
ligand_prmtop = app.AmberPrmtopFile('/Users/edwardkim/IGEM/gh_scoresolver/tests/pdb/ligandNOH.pdb')

# Create an OpenMM System from the ligand parameters
ligand_system = ligand_prmtop.createSystem()

# Serialize the System to an XML file
with open('ligand_forcefield.xml', 'w') as f:
    f.write(XmlSerializer.serialize(ligand_system))
