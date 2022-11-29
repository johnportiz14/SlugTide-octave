#!/bin/sh
'''
Add together the th0 and th90 degree strain results
to get areal strain.
- First, copy the corresponding azimuth files from directory above
'''
import numpy as np
from shutil import copyfile
import sys

# If this script called on its own (with no argument)
#---- (0th arg is script name)
if len(sys.argv)==1:
    print('Enter the name identifier used for the two files:')
    print('   (e.g., identifier for ``th90_welldw4`` is welldw4)')
    identifier = input()            # user enters the identifier
# If script called with commandline argument
else: identifier = sys.argv[1]


#  # Copy the azimuth files from directory above (working/) into this directory
#  copyfile('../th0_{}'.format(identifier),'th0_{}'.format(identifier))
#  print('successfully copied th0_{} into current directory.'.format(identifier))
#  copyfile('../th90_{}'.format(identifier),'th90_{}'.format(identifier))
#  print('successfully copied th90_{} into current directory.'.format(identifier))

# Read in 0° file first
d0 = np.genfromtxt('th0_{}'.format(identifier))
# Now read in 90° file
d90 = np.genfromtxt('th90_{}'.format(identifier))
# Add together to get areal strain, then save
areal_strain = d0+d90
#  np.savetxt('{}_tides.txt'.format(identifier), areal_strain)
np.savetxt('ertidAreal.txt', areal_strain)
np.savetxt('ertidAreal_{}.txt'.format(identifier), areal_strain)
