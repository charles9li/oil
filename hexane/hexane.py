from __future__ import division, print_function

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import gromacs
import mdtraj as md
import numpy as np

# load topology
top = gromacs.GromacsTopologyFile("hexane_L4nm.top")
top.box = (0., 0., 0., 3.9, 3.9, 3.9)*nanometers
top.box_vectors = np.diag([3.9, 3.9, 3.9])*nanometer

# create system
system = top.createSystem(nonbondedMethod=PME,
                          ewaldErrorTolerance=0.0001,
                          nonbondedCutoff=1.2*nanometer)
box_vectors = np.diag([3.9, 3.9, 3.9])*nanometer
system.setDefaultPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])

# integrator
temperature = 298.15*kelvin
friction = 1/picosecond
dt = 0.001*picosecond
integrator = LangevinIntegrator(temperature, friction, dt)

# barostat
pressure = 1*bar
barostat = MonteCarloBarostat(pressure, temperature)
system.addForce(barostat)

# create simulation
simulation = Simulation(top.topology, system, integrator)
t = md.load("hexane_L4nm.pdb")
simulation.context.setPositions(t.xyz[0])

# minimize energy
simulation.minimizeEnergy(maxIterations=1000)

# assign velocities
simulation.context.setVelocitiesToTemperature(temperature)

# equilibration
equilibration_steps = 20000
simulation.step(equilibration_steps)

# simulate
steps = 1000000
simulation.reporters.append(DCDReporter("hexane_298_01bar.dcd", 1000))
simulation.reporters.append(StateDataReporter("hexane_298_01bar.csv", 1000, time=True,
                                             potentialEnergy=True, kineticEnergy=True,
                                             volume=True, density=True, elapsedTime=True))
simulation.step(steps)
