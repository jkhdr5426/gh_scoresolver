from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from mdtraj.reporters import *
pdb = PDBFile("ahl/ahl_dock_luxr_1.pdb")
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml', 'ahl/gaff.xml', 'ahl/ahl-specific.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.deleteWater()
residues = modeller.addHydrogens(forcefield)
modeller.addSolvent(forcefield, padding=1.0*nanometer)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
print("Minimizing energy")
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
        potentialEnergy=True, temperature=True, volume=True, kineticEnergy=True, totalEnergy=True))
simulation.reporters.append(StateDataReporter("md_log.txt", 100, step=True,
        potentialEnergy=True, temperature=True, volume=True, kineticEnergy=True, totalEnergy=True))
simulation.reporters.append(NetCDFReporter('trajectory.nc', 100))
print("Running NVT")
simulation.step(10000)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)
print("Running NPT")
simulation.step(10000)
simulation.reporters[-1].close()


import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("md_log.txt", delimiter=',')

step = data[:,0]
potential_energy = data[:,1]
temperature = data[:,2]
volume = data[:,3]
kinetic_energy = data[:,4]
total_energy = data[:,5]

plt.plot(step, potential_energy, label="Potential Energy")
plt.plot(step, kinetic_energy, label="Kinetic Energy")
plt.plot(step, total_energy, label="Total Energy")
plt.xlabel("Step")
plt.ylabel("Energy (kJ/mol)")
plt.legend()
plt.show()

plt.plot(step, temperature)
plt.xlabel("Step")
plt.ylabel("Temperature (K)")
plt.show()
plt.plot(step, volume)
plt.xlabel("Step")
plt.ylabel("Volume (nm^3)")
plt.show()
