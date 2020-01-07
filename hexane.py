from __future__ import division, print_function

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import gromacs
import mdtraj

# load gro
gro = GromacsGroFile("hexane.gro")

# load topology
# top = GromacsTopFile("hexane.top", includeDir='~/gromacs/gromacs-5.1.1/share/top')
top = gromacs.GromacsTopologyFile("hexane.top")

# create system
system = top.createSystem()

# integrator
temperature = 300*kelvin
friction = 1/picosecond
dt = 0.002*picosecond
integrator = LangevinIntegrator(temperature, friction, dt)

# create simulation
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)

# minimize energy
simulation.minimizeEnergy(maxIterations=1000)

# equilibration
equilibration_steps = 1000
simulation.step(equilibration_steps)

# simulate
steps = 10000
simulation.reporters.append(PDBReporter("hexane.pdb", 100))
simulation._simulate(steps)
simulation.saveState("hexane_output.xml")
