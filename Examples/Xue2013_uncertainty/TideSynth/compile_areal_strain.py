'''
Add together the th0 and th90 degree strain results
to get areal strain.
- First, copy the corresponding azimuth files from directory above
'''
import numpy as np
from shutil import copyfile

print('Enter the name identifier used for the two files:')
print('   (e.g., identifier for ``th90_41-27`` is 41-27)')
identifier = input()            # user enters the identifier


# Copy the azimuth files from directory above into this directory
copyfile('../th0_{}'.format(identifier),'th0_{}'.format(identifier))
print('successfully copied th0_{} into current directory.'.format(identifier))
copyfile('../th90_{}'.format(identifier),'th90_{}'.format(identifier))
print('successfully copied th90_{} into current directory.'.format(identifier))

# Read in 0° file first
d0 = np.genfromtxt('th0_{}'.format(identifier))
# Now read in 90° file
d90 = np.genfromtxt('th90_{}'.format(identifier))
# Add together to get areal strain, then save
areal_strain = d0+d90
np.savetxt('{}_tides.txt'.format(identifier), areal_strain)
