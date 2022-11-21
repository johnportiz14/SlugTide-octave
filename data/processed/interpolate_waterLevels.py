#!/bin/sh
'''
A helper script to linearly interpolate a 2-column water level data set to a
certain time interval (e.g., to convert up-sample hourly data to 15-minute
intervals).
Only run this script after the datetimes have been formatted as serial dates in
Excel.
'''
import sys
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np

#------------------------------------------------------
#   EDIT USER PARAMETERS
#------------------------------------------------------
dt = 15.0               #[min] desired sample interval
num_header_rows = 0     # number of header rows to skip 
nan_string = 'NaT'      # missing value string (e.g., 'NaN')
#------------------------------------------------------


# If this script called on its own (with no argument)
#---- (0th arg is script name)
if len(sys.argv)==1:
    print('Enter the file name of the data you want to interpolate:')
    print('   (e.g.,  well_dw4.csv )')
    filename = input()            # user enters the filename
# If script called with commandline argument
else: filename = sys.argv[1]


if num_header_rows == 0: head = None
else: head = num_header_rows
df = pd.read_csv(filename,header=head,na_values=nan_string)

t = np.array(df.iloc[:,0])
wl = np.array(df.iloc[:,1])

#------------------------------------------------------
#  INTERPOLATE DATA
#------------------------------------------------------
print('Interpolating {} to {} min intervals...'.format(filename, dt))
# New times to interpolate WL to 
dti = dt/60./24. #convert dt to days
tnew = np.arange(t[0],t[-1], step=dti)
# Create interpolation object
f1 = interp1d(t, wl)
# Interpolate WL to new times
wlnew = f1(tnew)
print('  Done.')

#------------------------------------------------------
#  EXPORT DATA 
#------------------------------------------------------
new_filename = filename[:-4]+'_interp.csv'
print('Exporting to filename {}...'.format(new_filename))
twocols = np.column_stack( (tnew, wlnew) )
np.savetxt(new_filename, twocols, delimiter=',',fmt='%10.8f')
print('  Done.')


