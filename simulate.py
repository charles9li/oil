import argparse
import mdparse

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import gromacs
import mdtraj as md

import numpy as np


NONBONDED_METHODS = {'NoCutoff':            NonbondedForce.NoCutoff,
                     'CutoffNonPeriodic':   NonbondedForce.CutoffNonPeriodic,
                     'Ewald':               NonbondedForce.Ewald,
                     'PME':                 NonbondedForce.PME,
                     'LJPME':               NonbondedForce.LJPME}

BAROSTATS = [MonteCarloBarostat,
             MonteCarloAnisotropicBarostat,
             MonteCarloMembraneBarostat]


def _get_nonbonded_forces(system):
    """
    Returns a list of all nonbonded forces in the input system.

    :param system: OpenMM system
    :return: list of NonbondedForce instances
    """
    return [force for force in system.getForces() if isinstance(force, NonbondedForce)]


def _modify_barostat(system, ensemble_args):
    ensemble = ensemble_args.ensemble
    if ensemble in ['NVE', 'NVT']:
        _remove_barostat(system)
    if ensemble in ['NPT']:
        _add_barostat(system, ensemble_args)


def _add_barostat(system, ensemble_args):
    if not _has_barostat(system):
        temperature = ensemble_args.temperature*kelvin
        pressure = ensemble_args.pressure*bar
        barostat = MonteCarloBarostat(pressure, temperature)
        system.addForce(barostat)


def _remove_barostat(system):
    system_forces = system.getForces()
    for force in system_forces:
        if _is_barostat(force):
            i = system_forces.index(force)
            system_forces.pop(i)


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


def _create_simulation(topology, system, ensemble_args):
    integrator = _create_integrator(ensemble_args)
    simulation = Simulation(topology, system, integrator)
    return simulation


def _change_integrator(simulation, topology, system, ensemble_args):
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    box_vectors = state.getPeriodicBoxVectors()
    positions = state.getPositions()
    velocities = state.getVelocities()
    system.setDefaultPeriodicBoxVectors(*box_vectors)
    simulation = _create_simulation(topology, system, ensemble_args)
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)
    return simulation


def _create_integrator(ensemble_args):
    if ensemble_args.integrator == 'Langevin':
        temperature = ensemble_args.temperature*kelvin
        friction = ensemble_args.friction/picosecond
        time_step = ensemble_args.timestep*femtosecond
        return LangevinIntegrator(temperature, friction, time_step)
    if ensemble_args.integrator == 'Verlet':
        time_step = ensemble_args.timestep*femtosecond
        return VerletIntegrator(time_step)


def _set_positions(simulation, args):
    t = md.load(args.incoord)
    simulation.context.setPositions(t.xyz[0])


def _minimize_energy(simulation, ensemble_args):
    if ensemble_args.minimize:
        simulation.minimizeEnergy()


def _set_velocities_to_temperature(simulation, args):
    if args.gentemp > 0.0:
        simulation.context.setVelocitiesToTemperature(args.gentemp*kelvin)


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


# TODO: add ability to save states after running
def _run_simulation(simulation, ensemble_args):
    if not ensemble_args.average_volume and not ensemble_args.average_energy:
        simulation.step(ensemble_args.steps)
    else:
        volume_avg, energy_avg = _run_simulation_avg(simulation, ensemble_args)
        _save_state(simulation, ensemble_args)
        if ensemble_args.average_volume:
            simulation = _run_simulation_avg_vol(simulation, volume_avg)
            if ensemble_args.average_energy:
                _remove_barostat(simulation.system)
        if ensemble_args.average_energy:
            simulation = _run_simulation_avg_energy(simulation, energy_avg)
    _save_state(simulation, ensemble_args)
    return simulation


def _save_state(simulation, ensemble_args):
    if ensemble_args.savestate is not "":
        simulation.saveState(ensemble_args.savestate)


def _run_simulation_avg(simulation, ensemble_args):
    volume = []
    energy = []
    total_steps = ensemble_args.steps
    step_interval = 25  # TODO: need way to specify interval between taking volume and energy
    while total_steps > 0:
        simulation.step(step_interval)
        total_steps -= step_interval
        state = simulation.context.getState(getEnergy=True)
        volume.append(state.getPeriodicBoxVolume())
        energy.append(state.getKineticEnergy() + state.getPotentialEnergy())
    volume_avg = np.mean(volume)
    energy_avg = np.mean(energy)
    return volume_avg, energy_avg


def _run_simulation_avg_vol(simulation, volume_avg):
    step_interval = 25  # TODO: need way of specifying this
    while True:
        simulation.step(step_interval)
        volume_curr = simulation.context.getState().getPeriodicBoxVolume()
        if np.abs((volume_curr - volume_avg)/volume_avg) < 1e-3:
            break
    return simulation


def _run_simulation_avg_energy(simulation, energy_avg):
    step_interval = 25  # TODO: need way of specifying this
    while True:
        simulation.step(step_interval)
        state = simulation.context.getState(getEnergy=True)
        energy_curr = state.getKineticEnergy() + state.getPotentialEnergy()
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
    system = top.createSystem(nonbondedMethod=NoCutoff, ewaldErrorTolerance=args.ewald_error_tolerance,
                              nonbondedCutoff=args.nonbonded_cutoff*nanometer)
    for force in _get_nonbonded_forces(system):
        force.setUseDispersionCorrection(args.dispersion_correction)
        force.setNonbondedMethod(NONBONDED_METHODS[args.nonbonded_method])

    # create list of ensembles
    completed_ensembles = []

    # run simulation in each ensemble
    while args.ensembles:

        # remove ensemble from queue
        ensemble_args = args.ensembles.pop(0)

        # add or remove barostat
        _modify_barostat(system, ensemble_args)

        # initialize simulation
        if not completed_ensembles:
            simulation = _create_simulation(top.topology, system, ensemble_args)
            _set_positions(simulation, args)
            _set_velocities_to_temperature(simulation, args)
        else:
            simulation = _change_integrator(simulation, top.topology, system, ensemble_args)

        # minimize energy and equilibrate
        _minimize_energy(simulation, ensemble_args)
        _equilibrate(simulation, ensemble_args)

        # add and/or remove reporters
        _modify_reporters(simulation, ensemble_args)

        # run simulation
        simulation = _run_simulation(simulation, ensemble_args)

        # add ensemble to completed list
        completed_ensembles.append(ensemble_args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulation Options")
    parser.add_argument('parameter_file', default="params.in", help="Input file containing simulation parameters")
    command_line_args = parser.parse_args()

    # run
    main(command_line_args.parameter_file)
