# Written by Surilige Dhalenguite, Feb. 27th, 2019
# Physics College of Jilin University
# How to build a CIF file named "SrTiO3.cif"

import numpy as np

from ase import Atoms           # Add atoms, and you can add single atom with the help of "from ase import Atom"
from ase.build import bulk      # Build common bulk crystal
from ase.io import read, write  # Get I/O files with read/write
from ase.spacegroup import crystal, Spacegroup     # Import crystal parameters

import ase.io.cif               # Construct a CIF file

SrTiO3 = crystal(       # define SrTiO3 crystal
                 symbols=['Sr', 'Ti', 'O', 'O', 'O'],
                 basis=[(0.50000, 0.50000, 0.50000), # Sr1
                        (0.00000, 0.00000, 0.00000), # Ti1
                        (0.50000, 0.00000, 0.00000), # O1
                        (0.50000, 0.50000, 0.00000), # O2
                        (0.00000, 0.00000, 0.50000)], # O3
                 spacegroup='P1',
                 cellpar=[3.90500, 3.90500, 3.90500, 90, 90, 90])     # Crystal parameters

SrTiO3.write('SrTiO3.cif')
"""
    There may be several warnings after executing:
    
    UserWarning: scaled_positions 0 and 1 are equivalent 'are equivalent' % (kinds[ind], kind))
    UserWarning: scaled_positions 2 and 3 are equivalent 'are equivalent' % (kinds[ind], kind))
    
    But we'd still get the correct CIF file of FeSe. Fe1=atom(0), Fe2=atom(2), Se1=atom(3), Se2=atom(4), they are .
"""
