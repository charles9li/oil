from __future__ import division, print_function

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import parmed as pmd

topFile = 'Ethanol.top'
strucFile = 'Ethanol.gro'
top = pmd.load_file(topFile, xyz=strucFile)
