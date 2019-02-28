# Written by Surilige Dhalenguite, Feb. 27th, 2019
# Physics College of Jilin University
# How to build a CIF file named "FeSe_on_SrTiO3.cif" to observe the growth of FeSe crystal on SrTiO3 basis

import numpy as np

from ase import Atoms           # Add atoms, and you can add single atom with the help of "from ase import Atom"
from ase.build import bulk      # Build common bulk crystal
from ase.io import read, write  # Get I/O files with read/write
from ase.spacegroup import crystal, Spacegroup     # Import crystal parameters

from ase.build import add_adsorbate
from ase.calculators.emt import EMT

import ase.io.cif               # Construct a CIF file

FeSe_on_SrTiO3 = crystal(
                         symbols=['Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Sr', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'Ti', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'Fe', 'Fe', 'Fe', 'Fe', 'Fe', 'Fe', 'Fe', 'Fe', 'Se', 'Se', 'Se', 'Se', 'Se', 'Se', 'Se', 'Se'],
                         basis=[(0.25000, 0.25000, 0.07309), # Sr1
                                (0.75000, 0.25000, 0.07309), # Sr2
                                (0.25000, 0.75000, 0.07309), # Sr3
                                (0.75000, 0.75000, 0.07309), # Sr4
                                (0.25000, 0.25000, 0.21926), # Sr5
                                (0.75000, 0.25000, 0.21926), # Sr6
                                (0.25000, 0.75000, 0.21926), # Sr7
                                (0.75000, 0.75000, 0.21926), # Sr8
                                (0.25000, 0.25000, 0.36543), # Sr9
                                (0.75000, 0.25000, 0.36543), # Sr10
                                (0.25000, 0.75000, 0.36543), # Sr11
                                (0.75000, 0.75000, 0.36543), # Sr12
                                # There are 12 Sr atoms
                                (0.00000, 0.00000, 0.00000), # Ti1
                                (0.50000, 0.00000, 0.00000), # Ti2
                                (0.00000, 0.50000, 0.00000), # Ti3
                                (0.50000, 0.50000, 0.00000), # Ti4
                                (0.00000, 0.00000, 0.43852), # Ti5
                                (0.50000, 0.00000, 0.43852), # Ti6
                                (0.00000, 0.50000, 0.43852), # Ti7
                                (0.50000, 0.50000, 0.43852), # Ti8
                                (0.00000, 0.00000, 0.14617), # Ti9
                                (0.50000, 0.00000, 0.14617), # Ti10
                                (0.00000, 0.50000, 0.14617), # Ti11
                                (0.50000, 0.50000, 0.14617), # Ti12
                                (0.00000, 0.00000, 0.29235), # Ti13
                                (0.50000, 0.00000, 0.29235), # Ti14
                                (0.00000, 0.50000, 0.29235), # Ti15
                                (0.50000, 0.50000, 0.29235), # Ti16
                                # There are 16 Ti atoms
                                (0.25000, 0.00000, 0.00000), # O1
                                (0.00000, 0.25000, 0.00000), # O2
                                (0.25000, 0.00000, 0.43852), # O3
                                (0.00000, 0.25000, 0.43852), # O4
                                (0.00000, 0.00000, 0.07309), # O5
                                (0.75000, 0.00000, 0.00000), # O6
                                (0.50000, 0.25000, 0.00000), # O7
                                (0.75000, 0.00000, 0.43852), # O8
                                (0.50000, 0.25000, 0.43852), # O9
                                (0.50000, 0.00000, 0.07309), # O10
                                (0.25000, 0.50000, 0.00000), # O11
                                (0.00000, 0.75000, 0.00000), # O12
                                (0.25000, 0.50000, 0.43852), # O13
                                (0.00000, 0.75000, 0.43852), # O14
                                (0.00000, 0.50000, 0.07309), # O15
                                (0.75000, 0.50000, 0.00000), # O16
                                (0.50000, 0.75000, 0.00000), # O17
                                (0.75000, 0.50000, 0.43852), # O18
                                (0.50000, 0.75000, 0.43852), # O19
                                (0.50000, 0.50000, 0.07309), # O20
                                (0.25000, 0.00000, 0.14617), # O21
                                (0.00000, 0.25000, 0.14617), # O22
                                (0.00000, 0.00000, 0.21926), # O23
                                (0.75000, 0.00000, 0.14617), # O24
                                (0.50000, 0.25000, 0.14617), # O25
                                (0.50000, 0.00000, 0.21926), # O26
                                (0.25000, 0.50000, 0.14617), # O27
                                (0.00000, 0.75000, 0.14617), # O28
                                (0.00000, 0.50000, 0.21926), # O29
                                (0.75000, 0.50000, 0.14617), # O30
                                (0.50000, 0.75000, 0.14617), # O31
                                (0.50000, 0.50000, 0.21926), # O32
                                (0.25000, 0.00000, 0.29235), # O33
                                (0.00000, 0.25000, 0.29235), # O34
                                (0.00000, 0.00000, 0.36543), # O35
                                (0.75000, 0.00000, 0.29235), # O36
                                (0.25000, 0.50000, 0.29235), # O37
                                (0.50000, 0.00000, 0.36543), # O38
                                (0.25000, 0.50000, 0.29235), # O39
                                (0.00000, 0.75000, 0.29235), # O40
                                (0.00000, 0.50000, 0.36543), # O41
                                (0.75000, 0.50000, 0.29235), # O42
                                (0.50000, 0.75000, 0.29235), # O43
                                (0.50000, 0.00000, 0.36543), # O44
                                # There are 44 O atoms
                                (0.24139, 0.24139, 0.66457), # Fe1
                                (0.24139, 0.72417, 0.66457), # Fe2
                                (0.72417, 0.24139, 0.66457), # Fe3
                                (0.72417, 0.72417, 0.66457), # Fe4
                                (0.00000, 0.00000, 0.66457), # Fe5
                                (0.00000, 0.48278, 0.66457), # Fe6
                                (0.48278, 0.00000, 0.66457), # Fe7
                                (0.48278, 0.48278, 0.66457), # Fe8
                                # There are 8 Fe atoms
                                (0.24139, 0.00000, 0.71954), # Se1
                                (0.24139, 0.48278, 0.71954), # Se2
                                (0.72417, 0.00000, 0.71954), # Se3
                                (0.72417, 0.48278, 0.71954), # Se4
                                (0.00000, 0.24139, 0.60960), # Se5
                                (0.00000, 0.72417, 0.60960), # Se6
                                (0.48278, 0.24139, 0.60960), # Se7
                                (0.48278, 0.72417, 0.60960)],# Se8
                                # There are 8 Se atoms
                         spacegroup='P1',  # nor replaced by corresponding figure
                         cellpar=[7.81000, 7.81000, 26.71500, 90, 90, 90])     # Crystal parameters

FeSe_on_SrTiO3.write('FeSe_on_SrTiO3.cif')

"""
    There may be several warnings after executing:
    
    UserWarning: scaled_positions 0 and 1 are equivalent 'are equivalent' % (kinds[ind], kind))
    UserWarning: scaled_positions 2 and 3 are equivalent 'are equivalent' % (kinds[ind], kind))
    
    But we'd still get the correct CIF file of FeSe. Fe1=atom(0), Fe2=atom(2), Se1=atom(3), Se2=atom(4), they are .
"""
