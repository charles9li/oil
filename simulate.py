import argparse
import mdparse

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from parmed import gromacs


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


def _create_integrator(ensemble_args):
    if ensemble_args.integrator == 'Langevin':
        temperature = ensemble_args.temperature*kelvin
        friction = ensemble_args.friction/picosecond
        time_step = ensemble_args.timestep*femtosecond
        return LangevinIntegrator(temperature, friction, time_step)
    if ensemble_args.integrator == 'Verlet':
        time_step = ensemble_args.timestep*femtosecond
        return VerletIntegrator(time_step)


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
    system_forces = system.getForces
    for force in system_forces:
        if _is_barostat(force):
            i = system_forces.index(force)
            system_forces.pop(i)
            return


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


def main(parameter_file="params.in"):
    """

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
        ensemble_args = args.ensembles.pop(0)
        integrator = _create_integrator(ensemble_args)
        _modify_barostat(system, ensemble_args)
        if not completed_ensembles:
            pass
        else:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulation Options")
    parser.add_argument('parameter_file', default="params.in", help="Input file containing simulation parameters")
    command_line_args = parser.parse_args()

    # run
    main(command_line_args.parameter_file)
