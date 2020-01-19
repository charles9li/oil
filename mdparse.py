from collections import OrderedDict
from ast import literal_eval


class _Options(object):

    def __getattr__(self, attr):
        return self.activeOptions[attr]

    def __init__(self):
        self.userOptions = OrderedDict()
        self.activeOptions = OrderedDict()
        self.documentation = OrderedDict()

    def set_active(self, key, default, value_type, doc="No documentation available.", allowed=None):
        """ Set one option.  The arguments are:
        key     : The name of the option.
        default : The default value.
        typ     : The type of the value.
        doc     : The documentation string.
        allowed : An optional list of allowed values.
        depend  : A condition that must be True for the option to be activated.
        clash   : A condition that must be False for the option to be activated.
        msg     : A warning that is printed out if the option is not activated.
        """
        # find value of key
        if key in self.userOptions:
            try:
                val = literal_eval(self.userOptions[key])
            except ValueError:
                val = str(self.userOptions[key])
        else:
            val = default

        # check if value is allowed
        if type(allowed) is list:
            if val not in allowed:
                raise Exception("'{0}' is not an allowed for for the '{1}' parameter.".format(val, key))

        # check if value is of correct type
        if val is not None and type(val) is not value_type:
            raise Exception("'{0}' cannot be of '{1}' type.".format(val, type(val).__name__))

        # set value to key in activeOptions
        self.activeOptions[key] = val

        # adds documentation
        self.documentation[key] = doc


class SimulationOptions(_Options):

    def __init__(self, input_file):
        super(SimulationOptions, self).__init__()

        self.ensembles = []

        # read input file and gather parameter/value pairs
        option_parser = self
        if input_file is not None:
            with open(input_file) as f:
                for line in f:
                    line = line.strip().partition('#')[0]
                    line = line.rstrip()
                    s = line.split()
                    if len(s) > 0:
                        if len(s) > 2:
                            raise Exception("Only 1 argument can be specified. {:d} were given.".format(len(s) - 1))
                        key = s[0]
                        if key == 'ensemble':
                            option_parser = _EnsembleOptions()
                            self.ensembles.append(option_parser)
                        try:
                            val = s[1]
                        except IndexError:
                            val = None
                        option_parser.userOptions[key] = val

        # set parameters for input files
        self.set_active('topfile', 'system.top', str,
                        doc="Gromacs system.top file")
        self.set_active('grofile', 'box.gro', str,
                        doc="Gromacs .gro file, we just use for box")
        self.set_active('incoord', 'in.pdb', str,
                        doc="input file, must be .pdb or .xml")

        # set parameters for force field
        self.set_active('nonbonded_method', 'PME', str,
                        doc="Set the method for nonbonded interactions.",
                        allowed=["NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "PME", "LJPME"])
        self.set_active('nonbonded_cutoff', 0.9, float,
                        doc="Nonbonded cutoff distance in nanometers.")
        self.set_active('dispersion_correction', True, bool,
                        doc="Isotropic long-range dispersion correction for periodic systems.")
        self.set_active('ewald_error_tolerance', 0.0005, float,
                        doc="Error tolerance for Ewald, PME, LJPME methods. Don't go below 5e-5 for PME unless running "
                            "in double precision.")
        self.set_active('dispersion_correction', True, bool,
                        doc="Isotropic long-range dispersion correction for periodic systems.")

        # set simulation conditions
        self.set_active('temperature', 0.0, float,
                        doc="Simulation temperature for Langevin integrator or Andersen thermostat.")
        self.set_active('gentemp', self.temperature, float,
                        doc="Specify temperature for generating velocities")
        self.set_active('pressure', 0.0, float,
                        doc="Simulation pressure; set a positive number to activate.")
        self.set_active('collision_rate', 1.0, float,
                        "Collision rate of system is 1/ps, used in some integrators.")

        # set parameters for each ensemble
        for ensemble_args in self.ensembles:
            ensemble_args.set_active_all()
            if ensemble_args.ensemble in ['NVT', 'NPT']:
                ensemble_args.userOptions['temperature'] = str(self.temperature)
                ensemble_args.set_active('temperature', self.temperature, float,
                                         doc="Temperature of simulation")
            if ensemble_args.ensemble in ['NPT']:
                ensemble_args.userOptions['pressure'] = str(self.pressure)
                ensemble_args.set_active('pressure', self.pressure, float,
                                         doc="Pressure of simulation")
            if ensemble_args.integrator in ['Langevin']:
                ensemble_args.userOptions['friction'] = str(self.collision_rate)
                ensemble_args.set_active('friction', 1.0, float,
                                         doc="Friction used in some integrators.")


class _EnsembleOptions(_Options):

    ALLOWED_INTEGRATORS = {'NVE': ['Verlet'],
                           'NVT': ['Langevin'],
                           'NPT': ['Langevin']}

    def __init__(self):
        super(_EnsembleOptions, self).__init__()

    def set_active_all(self):

        # set simulation parameters
        self.set_active('ensemble', 'NVE', str,
                        doc="Ensemble",
                        allowed=['NVE', 'NVT', 'NPT'])
        self.set_active('integrator', 'Verlet', str,
                        doc="Integrator to use for simulation",
                        allowed=self.ALLOWED_INTEGRATORS[self.ensemble])
        self.set_active('timestep', 1.0, float,
                        doc="Timestep in femtoseconds.")
        self.set_active('minimize', False, bool,
                        doc="Whether or not to minimize energy")
        self.set_active('equilibrate', 0, int,
                        doc="Number of steps reserved for equilibration")
        self.set_active('steps', 1000000, int,
                        doc="Number of steps for simulation")
        self.set_active('outstate', "output.csv", str,
                        doc="state output file name")
        self.set_active('outdcd', "output.dcd", str,
                        doc="DCD output file")
        self.set_active('state_report_interval', 0, int,
                        doc="Number of steps between every state report")
        self.set_active('dcd_report_interval', 0, int,
                        doc="Number of steps between every DCD report")
        self.set_active('savestate', "", str,
                        'xml file that contains state information')
        self.set_active('average_volume', False, bool,
                        doc="Whether or not to run simulation so that it averages volume.")
        self.set_active('average_energy', False, bool,
                        doc="Whether or not to run simulation so that it averages energy.")


if __name__ == '__main__':
    args = SimulationOptions("params.in")
