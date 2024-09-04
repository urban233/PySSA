"""This is a working OpenMM simulation script."""
from sys import stdout

from openmm.app import *
from openmm import *
from openmm.unit import *

# 1.loading initial coordinates
pdb = PDBFile(r"C:\Users\student\sciebo\ProteinStructurePrediction\Protein Modelling\Proteins\pdb\AlphaFold 3\BMP2-8M_alphafold3.pdb")

# 2.choosing a forcefield parameters
ff = ForceField('amber10.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(ff)

system = ff.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic)

# 3. Choose parameters of the experiment: temperature, pressure, box size, solvation, boundary conditions, etc
temperature = 300*kelvin
frictionCoeff = 1/picosecond
time_step = 0.002*picoseconds
total_steps = 400*picoseconds / time_step

# 4. Choose an algorithm (integrator)
integrator = LangevinIntegrator(temperature, frictionCoeff, time_step)

# 5. Run simulation, saving coordinates time to time:
# 5a. Create a simulation object
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
# 5b. Minimize energy
simulation.minimizeEnergy()
# 5c. Save coordinates to dcd file and energies to a standard output console:
simulation.reporters.append(
  DCDReporter(file='example.dcd', reportInterval=1000)
)
simulation.reporters.append(
  StateDataReporter(
    stdout, 5000, step=True, potentialEnergy=True, temperature=True, progress=True, totalSteps = total_steps)
)
# 5d. Run!
simulation.step(total_steps)

print('Saving...')
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))
print('Done')
