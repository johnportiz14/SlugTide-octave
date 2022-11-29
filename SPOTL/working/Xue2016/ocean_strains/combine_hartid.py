'''
Add together the hart0 and hart90 degree harmonic components.
'''
import numpy as np

# Read in 0° file first
#  d0 = np.genfromtxt('th0_{}'.format(identifier))
d0 = np.genfromtxt('hart0')
# Now read in 90° file
#  d90 = np.genfromtxt('th90_{}'.format(identifier))
d90 = np.genfromtxt('hart90')
# Add together to get combined areal harmonics, then save
areal_strain = d0+d90
np.savetxt('hartidAreal.txt', areal_strain)

