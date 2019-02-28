# Written by Surilige Dhalenguite, Feb. 27th, 2019
# Physics College of Jilin University
# How to build a CIF file named "FeSe.cif"

import numpy as np

from ase import Atoms           # Add atoms, and you can add single atom with the help of "from ase import Atom"
from ase.build import bulk      # Build common bulk crystal
from ase.io import read, write  # Get I/O files with read/write
from ase.spacegroup import crystal, Spacegroup     # Import crystal parameters
# a = 3.77050, b = 3.77050, c = 5.50800 are crystal parameters for FeSe bulk

import ase.io.cif               # Construct a CIF file
# Atoms = ase.io.cif.read_cif("FeSe.cif")

FeSe = crystal(     # define FeSe crystal
               symbols=['Fe', 'Fe', 'Se', 'Se'],
               basis=[(0.50000, 0.50000, 0.50000), # Fe1
                      (0.00000, 0.00000, 0.50000), # Fe2
                      (0.50000, 0.00000, 0.76660), # Se1
                      (0.00000, 0.50000, 0.23340)], # Se2
               spacegroup='P4/nmm',  # nor replaced by corresponding figure
               cellpar=[3.77050, 3.77050, 5.50800, 90, 90, 90])     # Crystal parameters

FeSe.write('FeSe.cif')
"""
    There may be several warnings after executing:
    
    UserWarning: scaled_positions 0 and 1 are equivalent 'are equivalent' % (kinds[ind], kind))
    UserWarning: scaled_positions 2 and 3 are equivalent 'are equivalent' % (kinds[ind], kind))
    
    But we'd still get the correct CIF file of FeSe. Fe1=atom(0), Fe2=atom(2), Se1=atom(3), Se2=atom(4), they are .
"""
