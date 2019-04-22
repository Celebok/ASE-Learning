# Written by Surilige Dhalenguite, Apr. 17th, 2019
# Physics College of Jilin University
# How to build a CIF file named "FeSe.cif"

import numpy as np

from ase import Atoms           # Add atoms, and you can add single atom with the help of "from ase import Atom"
from ase.build import bulk      # Build common bulk crystal
from ase.io import read, write  # Get I/O files with read/write
from ase.spacegroup import crystal, Spacegroup     # Import crystal parameters

import ase.io.cif               # Construct a CIF file

La_Fmmm = crystal(     # define La crystal
               symbols=['La', 'La', 'La', 'La'],
	       basis=[(0.00000, 0.50000, 0.00000), # La1
                      (0.00000, 0.00000, 0.50000), # La2
		      (0.50000, 0.50000, 0.50000), # La3
		      (0.50000, 0.00000, 0.00000)], # La4
               spacegroup='Fmmm',  # nor replaced by corresponding figure
               cellpar=[3.94790, 4.07160, 3.99010, 90, 90, 90])     # Crystal parameters
La_Fmmm.write('La_Fmmm.cif')

"""
 There may be several warnings after executing:
 
 UserWarning: scaled_positions 0 and 1 are equivalent 'are equivalent' % (kinds[ind], kind))
 UserWarning: scaled_positions 0 and 2 are equivalent 'are equivalent' % (kinds[ind], kind))
 UserWarning: scaled_positions 0 and 3 are equivalent 'are equivalent' % (kinds[ind], kind))    

 But there is nothing.
"""
