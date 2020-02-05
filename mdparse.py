from collections import OrderedDict
from ast import literal_eval
import argparse


class _ParameterInfo(object):

    def __init__(self, default_val, val_type, doc=None, allowed=None):
        self.default = default_val
        self.type = val_type
        self.documentation = doc
        self.allowed = allowed


_ALLOWED_NONBONDED_METHODS = ["NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "PME", "LJPME"]

_SIMULATION_PARAMETER_INFO = {'topfile':                _ParameterInfo('system.top', str),
                              'grofile':                _ParameterInfo('box.gro', str),
                              'incoord':                _ParameterInfo('in.pdb', str),
                              'nonbonded_method':       _ParameterInfo('PME', str, allowed=_ALLOWED_NONBONDED_METHODS),
                              'nonbonded_cutoff':       _ParameterInfo(0.9, float),
                              'dispersion_correction':  _ParameterInfo(True, bool),
                              'ewald_error_tolerance':  _ParameterInfo(0.005, float),
                              'temperature':            _ParameterInfo(0.0, float),
                              'gentemp':                _ParameterInfo(0.0, float),
                              'pressure':               _ParameterInfo(0.0, float),
                              'collision_rate':         _ParameterInfo(1.0, float)}

_ALLOWED_INTEGRATORS = {'NVE': ['Verlet'],
                        'NVT': ['Langevin'],
                        'NPT': ['Langevin']}
_ALLOWED_ENSEMBLES = list(_ALLOWED_INTEGRATORS.keys())

_ENSEMBLE_PARAMETER_INFO = {'ensemble':                 _ParameterInfo('NVE', str, allowed=_ALLOWED_ENSEMBLES),
                            'integrator':               _ParameterInfo('Verlet', str),
                            'timestep':                 _ParameterInfo(1.0, float),
                            'minimize':                 _ParameterInfo(False, bool),
                            'equilibrate':              _ParameterInfo(0, int),
                            'steps':                    _ParameterInfo(1000000, int),
                            'outstate':                 _ParameterInfo("output.csv", str),
                            'outdcd':                   _ParameterInfo("output.dcd", str),
                            'state_report_interval':    _ParameterInfo(0, int),
                            'dcd_report_interval':      _ParameterInfo(0, int),
                            'savestate':                _ParameterInfo("output.xml", str),
                            'average_volume':           _ParameterInfo(False, True),
                            'average_energy':           _ParameterInfo(False, bool)}


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
                raise Exception("'{0}' is not an allowed for for the '{1}' parameter.".format(val, key))

        # check if value is of correct type
        if val is not None and type(val) is not value_type:
            raise Exception("'{0}' cannot be of '{1}' type.".format(val, type(val).__name__))

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
        for key, val in _SIMULATION_PARAMETER_INFO:
            self.set_active(key, val.default, val.type, doc=val.doc, allowed=val.allowed)

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
        for key, val in _ENSEMBLE_PARAMETER_INFO:
            if key == 'integrator':
                self.set_active(key, val.default, val.type, doc=val.doc, allowed=_ALLOWED_INTEGRATORS[self.ensemble])
            else:
                self.set_active(key, val.default, val.type, doc=val.doc, allowed=val.allowed)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.parse_args()
