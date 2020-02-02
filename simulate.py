"""
The simulate module is used to run an OpenMM simulation. The mdparse module is called to parse user
inputs from an input file and determine simulation conditions, force methods, and integrators.

@author Charles Li <charlesli@ucsb.edu>
"""

# imports for command line arguments
import argparse
import mdparse

# OpenMM imports
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import gromacs

# MDTraj import
import mdtraj as md

# Numpy import
import numpy as np

# Custom statistical analysis module
from openmm_stats import *

# Dictionary of nonbonded methods
NONBONDED_METHODS = {'NoCutoff':            NonbondedForce.NoCutoff,
                     'CutoffNonPeriodic':   NonbondedForce.CutoffNonPeriodic,
                     'Ewald':               NonbondedForce.Ewald,
                     'PME':                 NonbondedForce.PME,
                     'LJPME':               NonbondedForce.LJPME}

# List of barostats
BAROSTATS = [MonteCarloBarostat,
             MonteCarloAnisotropicBarostat,
             MonteCarloMembraneBarostat]


def _create_system(top, args):
    system = top.createSystem(nonbondedMethod=NoCutoff, ewaldErrorTolerance=args.ewald_error_tolerance,
                              nonbondedCutoff=args.nonbonded_cutoff*nanometer)
    for force in _get_nonbonded_forces(system):
        force.setUseDispersionCorrection(args.dispersion_correction)
        force.setNonbondedMethod(NONBONDED_METHODS[args.nonbonded_method])
    return system


def _get_nonbonded_forces(system):
    """
    Returns a list of all nonbonded forces in the input system.

    :param system: OpenMM system
    :return: list of NonbondedForce instances
    """
    return [force for force in system.getForces() if isinstance(force, NonbondedForce)]


def _modify_barostat(system, ensemble_args):
    """
    Modifies the barostat depending on the ensemble in which the current simulation is running in.
    Adds a barostat if running NPT and removes it if running NVE or NVT.

    :param system: OpenMM system
    :param ensemble_args: ensemble arguments
    """
    ensemble = ensemble_args.ensemble
    if ensemble in ['NVE', 'NVT']:
        _remove_barostat(system)
    if ensemble in ['NPT']:
        _add_barostat(system, ensemble_args)


def _add_barostat(system, ensemble_args):
    if not _has_barostat(system):
        temperature = ensemble_args.temperature*kelvin
        pressure = ensemble_args.pressure*bar
        # TODO: add functionality for other barostats
        # TODO: add option to specify barostat interval
        barostat = MonteCarloBarostat(pressure, temperature)
        system.addForce(barostat)


def _remove_barostat(system):
    system_forces = system.getForces()
    for force in system_forces:
        if _is_barostat(force):
            i = system_forces.index(force)
            system.removeForce(i)


def _has_barostat(system):
    for force in system.getForces():
        if _is_barostat(force):
            return True
    return False


def _is_barostat(force):
    for barostat in BAROSTATS:
        if isinstance(force, barostat):
            return True
    return False


def _create_simulation(simulation, topology, system, args, ensemble_args):
    # retrieve positions, velocities, and box vectors
    if simulation is None:
        positions = _get_positions_from_incoord(args)
        velocities = None
        box_vectors = None
    else:
        state = simulation.context.getState(getPositions=True, getVelocities=True)
        positions = state.getPositions()
        velocities = state.getVelocities()
        box_vectors = state.getPeriodicBoxVectors()

    # create integrator
    integrator = _create_integrator(ensemble_args)

    # create simulation
    simulation = Simulation(topology, system, integrator)

    # set positions
    simulation.context.setPositions(positions)

    # set velocities
    if velocities is None:
        _set_velocities_to_temperature(simulation, args)
    else:
        simulation.context.setVelocities(velocities)

    # set box vectors
    if box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*box_vectors)

    return simulation


def _get_positions_from_incoord(args):
    t = md.load(args.incoord)
    return t.xyz[0]


def _create_integrator(ensemble_args):
    if ensemble_args.integrator == 'Langevin':
        temperature = ensemble_args.temperature*kelvin
        friction = ensemble_args.friction/picosecond
        time_step = ensemble_args.timestep*femtosecond
        return LangevinIntegrator(temperature, friction, time_step)
    if ensemble_args.integrator == 'Verlet':
        time_step = ensemble_args.timestep*femtosecond
        return VerletIntegrator(time_step)


def _set_velocities_to_temperature(simulation, args):
    if args.gentemp > 0.0:
        simulation.context.setVelocitiesToTemperature(args.gentemp*kelvin)


def _minimize_energy(simulation, ensemble_args):
    if ensemble_args.minimize:
        simulation.minimizeEnergy()


def _equilibrate(simulation, ensemble_args):
    equilibration_steps = ensemble_args.equilibrate
    if equilibration_steps > 0:
        simulation.step(equilibration_steps)


def _modify_reporters(simulation, ensemble_args):

    # remove current reporters
    while simulation.reporters:
        simulation.reporters.pop()

    # add reporters
    _add_state_reporter(simulation, ensemble_args)
    _add_dcd_reporter(simulation, ensemble_args)


def _add_state_reporter(simulation, ensemble_args):
    report_interval = ensemble_args.state_report_interval
    if report_interval > 0:
        filename = ensemble_args.outstate
        state_data_reporter = StateDataReporter(filename, report_interval,
                                                step=True, time=True,
                                                potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                                temperature=True, volume=True, density=True,
                                                speed=True, elapsedTime=True)
        simulation.reporters.append(state_data_reporter)


def _add_dcd_reporter(simulation, ensemble_args):
    report_interval = ensemble_args.dcd_report_interval
    if report_interval > 0:
        filename = ensemble_args.outdcd
        dcd_reporter = DCDReporter(filename, report_interval, enforcePeriodicBox=False)
        simulation.reporters.append(dcd_reporter)


def _run_simulation(simulation, args, ensemble_args):
    # Check that StateDataReporter is appended to simulation if averaging volume and/or energy
    if ensemble_args.average_volume or ensemble_args.average_energy:
        if ensemble_args.state_report_interval == 0:
            raise ValueError("To average volume and/or energy, there must be a StateDataReporter "
                             "appended to the simulation.")

    # Step simulation
    simulation.step(ensemble_args.steps)

    # Average volume and/or energy if specified
    if ensemble_args.average_volume or ensemble_args.average_energy:
        volume_avg = compute_stats_from_file_by_column_name(ensemble_args.outstate, 'Box Volume (nm^3)')[2]
        energy_avg = compute_stats_from_file_by_column_name(ensemble_args.outstate, 'Total Energy (kJ/mole)')[2]
        _save_state(simulation, ensemble_args)
        if ensemble_args.average_volume:
            simulation = _run_simulation_avg_vol(simulation, volume_avg)
            if ensemble_args.average_energy:
                simulation = _remove_barostat_from_simulation(simulation, args, ensemble_args)
        if ensemble_args.average_energy:
            simulation = _run_simulation_avg_energy(simulation, energy_avg)
    _save_state(simulation, ensemble_args)
    return simulation


def _remove_barostat_from_simulation(simulation, args, ensemble_args):
    system = simulation.system
    _remove_barostat(system)
    topology = simulation.topology
    simulation = _create_simulation(simulation, topology, system, args, ensemble_args)
    return simulation


def _save_state(simulation, ensemble_args):
    if ensemble_args.savestate is not "":
        simulation.saveState(ensemble_args.savestate)


def _run_simulation_avg_vol(simulation, volume_avg):
    step_interval = 25  # TODO: need way of specifying this
    while True:
        simulation.step(step_interval)
        volume_curr = simulation.context.getState().getPeriodicBoxVolume()
        # TODO: have way to specify error
        if np.abs((volume_curr - volume_avg)/volume_avg) < 1e-3:
            break
    return simulation


def _run_simulation_avg_energy(simulation, energy_avg):
    step_interval = 25  # TODO: need way of specifying this
    while True:
        simulation.step(step_interval)
        state = simulation.context.getState(getEnergy=True)
        energy_curr = state.getKineticEnergy() + state.getPotentialEnergy()
        # TODO: have way to specify error
        if np.abs((energy_curr - energy_avg)/energy_avg) < 1e-3:
            break
    return simulation


def _ensemble_name_list(args):
    return [ensemble_args.ensemble for ensemble_args in args.ensembles]


def main(parameter_file="params.in"):
    """
    Creates an OpenMM system and runs simulations with the appropriate integrator using
    specifications in an input file.

    :param parameter_file: input file containing parameters
    :return: None
    """
    # parse user input
    args = mdparse.SimulationOptions(parameter_file)

    # load topology
    top = gromacs.GromacsTopologyFile(args.topfile)
    gro = gromacs.GromacsGroFile.parse(args.grofile)
    top.box = gro.box

    # create system
    system = _create_system(top, args)

    # initialize simulation
    simulation = None

    # run simulation in each ensemble
    while args.ensembles:

        # remove ensemble from queue
        ensemble_args = args.ensembles.pop(0)

        # add barostat if in NPT
        _modify_barostat(system, ensemble_args)

        # initialize simulation
        simulation = _create_simulation(simulation, top.topology, system, args, ensemble_args)

        # minimize energy and equilibrate
        _minimize_energy(simulation, ensemble_args)
        _equilibrate(simulation, ensemble_args)

        # add and/or remove reporters
        _modify_reporters(simulation, ensemble_args)

        # run simulation
        simulation = _run_simulation(simulation, args, ensemble_args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulation Options")
    parser.add_argument('parameter_file', default="params.in", help="Input file containing simulation parameters")
    command_line_args = parser.parse_args()

    # run
    main(command_line_args.parameter_file)
