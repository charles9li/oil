"""
The mdparse module is used to parse user input from an input file. The information is stored in a
SimulationOptions object and passed to simulate.py to run an OpenMM simulation.

@author Charles Li <charlesli@ucsb.edu>
"""
from collections import OrderedDict
from ast import literal_eval


class _ParameterInfo(object):
    """
    This class contains the following information for each parameter:
        1. default value
        2. type
        3. documentation
        4. allowed values
    """

    def __init__(self, default_val, val_type, doc="No documentation available.", allowed=None):
        """
        Create a _ParameterInfo instance.

        Parameters
        ----------
        default_val : val_type
            default value of parameter
        val_type : type
            parameter type
        doc : str
            parameter documentation
        allowed : NoneType or list
            allowed values for parameter, if applicable
        """
        self.default = default_val
        self.type = val_type
        self.doc = doc
        self.allowed = allowed


# allowed values for nonbonded methods, integrators, and ensembles
_ALLOWED_NONBONDED_METHODS = ["NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "PME", "LJPME"]
_ALLOWED_INTEGRATORS = {'NVE': ['Verlet'],
                        'NVT': ['Langevin', 'NoseHooverChainVelocityVerlet'],
                        'NPT': ['Langevin', 'NoseHooverChainVelocityVerlet']}
_ALLOWED_ENSEMBLES = list(_ALLOWED_INTEGRATORS.keys())

# parameter information for the main simulation
_SIMULATION_PARAMETER_INFO = OrderedDict()
_SIMULATION_PARAMETER_INFO['topfile'] =                 _ParameterInfo('system.top', str)
_SIMULATION_PARAMETER_INFO['grofile'] =                 _ParameterInfo('box.gro', str)
_SIMULATION_PARAMETER_INFO['incoord'] =                 _ParameterInfo('in.pdb', str)
_SIMULATION_PARAMETER_INFO['nonbonded_method'] =        _ParameterInfo('PME', str, allowed=_ALLOWED_NONBONDED_METHODS)
_SIMULATION_PARAMETER_INFO['nonbonded_cutoff'] =        _ParameterInfo(0.9, float)
_SIMULATION_PARAMETER_INFO['dispersion_correction'] =   _ParameterInfo(True, bool)
_SIMULATION_PARAMETER_INFO['ewald_error_tolerance'] =   _ParameterInfo(0.005, float)
_SIMULATION_PARAMETER_INFO['temperature'] =             _ParameterInfo(0.0, float)
_SIMULATION_PARAMETER_INFO['gentemp'] =                 _ParameterInfo(0.0, float)
_SIMULATION_PARAMETER_INFO['pressure'] =                _ParameterInfo(0.0, float)
_SIMULATION_PARAMETER_INFO['collision_rate'] =          _ParameterInfo(1.0, float)

# parameter information for each individual ensemble
_ENSEMBLE_PARAMETER_INFO = OrderedDict()
_ENSEMBLE_PARAMETER_INFO['ensemble'] =              _ParameterInfo('NVE', str, allowed=_ALLOWED_ENSEMBLES)
_ENSEMBLE_PARAMETER_INFO['integrator'] =            _ParameterInfo('Verlet', str)
_ENSEMBLE_PARAMETER_INFO['timestep'] =              _ParameterInfo(1.0, float)
_ENSEMBLE_PARAMETER_INFO['minimize'] =              _ParameterInfo(False, bool)
_ENSEMBLE_PARAMETER_INFO['equilibrate'] =           _ParameterInfo(0, int)
_ENSEMBLE_PARAMETER_INFO['steps'] =                 _ParameterInfo(1000000, int)
_ENSEMBLE_PARAMETER_INFO['outstate'] =              _ParameterInfo("output.csv", str)
_ENSEMBLE_PARAMETER_INFO['outdcd'] =                _ParameterInfo("output.dcd", str)
_ENSEMBLE_PARAMETER_INFO['state_report_interval'] = _ParameterInfo(0, int)
_ENSEMBLE_PARAMETER_INFO['dcd_report_interval'] =   _ParameterInfo(0, int)
_ENSEMBLE_PARAMETER_INFO['savestate'] =             _ParameterInfo("output.xml", str)
_ENSEMBLE_PARAMETER_INFO['average_volume'] =        _ParameterInfo(False, bool)
_ENSEMBLE_PARAMETER_INFO['average_energy'] =        _ParameterInfo(False, bool)


class _Options(object):
    """
    The Options class has methods to set values for parameters. Parent class of SimulationOptions
    and _EnsembleOptions.
    """

    def __getattr__(self, attr):
        return self.activeOptions[attr]

    def __init__(self):
        self.userOptions = OrderedDict()
        self.activeOptions = OrderedDict()
        self.documentation = OrderedDict()

    def set_active(self, key, default, value_type, doc="No documentation available.", allowed=None):
        """ Set one option.  The arguments are:

        Parameters
        ----------
        key : str
            parameter name
        default : value_type
            default value of parameter
        value_type : type
            parameter type
        doc : str
            parameter documentation
        allowed : NoneType or list
            allowed values for parameter, if applicable
        """
        # find value of key
        is_default = False
        if key in self.userOptions:
            try:
                val = literal_eval(self.userOptions[key])
            except ValueError:
                val = str(self.userOptions[key])
        else:
            val = default
            is_default = True

        # check if value is allowed
        if type(allowed) is list:
            if val not in allowed:
                raise Exception("'{}' is not an allowed for for the '{}' parameter.".format(val, key))

        # check if value is of correct type
        if val is not None and type(val) is not value_type:
            raise Exception("For the '{}' parameter, '{}' cannot be type '{}'.".format(key, val, value_type.__name__))

        # set value to key in activeOptions
        self.activeOptions[key] = val
        if is_default:
            self._print_option_default(key, val)
        else:
            self._print_option_set(key, val)

        # adds documentation
        self.documentation[key] = doc

    @staticmethod
    def _print_option_set(key, val):
        print("'{0}' parameter set to value of '{1}'.".format(key, val))

    @staticmethod
    def _print_option_default(key, val):
        print(("No argument specified for the '{0}' parameter. "
               "Setting '{0}' parameter to default value of '{1}'.").format(key, val))


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

        # print header
        print("\nParsing parameters from input file.")
        print("===================================")

        # set simulation parameters
        for parameter in _SIMULATION_PARAMETER_INFO:
            info = _SIMULATION_PARAMETER_INFO[parameter]
            self.set_active(parameter, info.default, info.type, doc=info.doc, allowed=info.allowed)

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

        # print header
        print("\nParsing ensemble arguments.")
        print("===========================")

        # set ensemble options
        for parameter in _ENSEMBLE_PARAMETER_INFO:
            info = _ENSEMBLE_PARAMETER_INFO[parameter]
            if parameter == 'integrator':
                allowed = _ALLOWED_INTEGRATORS[self.ensemble]
            else:
                allowed = info.allowed
            self.set_active(parameter, info.default, info.type, doc=info.doc, allowed=allowed)
