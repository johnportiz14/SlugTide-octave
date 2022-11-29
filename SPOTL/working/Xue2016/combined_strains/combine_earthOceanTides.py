'''
Add together the areal strains from both Earth and Ocean tides.
'''
import numpy as np

# Read in Ocean tides first 
do = np.genfromtxt('hartidAreal.txt')
# Now read in Earth tides 
de = np.genfromtxt('ertidAreal.txt')
# Add together to get combined Earth+Ocean strains 
combined_areal_strain = do+de
np.savetxt('combinedArealStrains.txt', combined_areal_strain)

