#!/bin/sh
'''
Convert local times (which may be in Daylight Savings or Standard Time to UTC
time.
'''
import sys,os
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)
import pandas as pd

#------------------------------------------------------
#   EDIT USER PARAMETERS
#------------------------------------------------------
#  data_dir = '/project/gas_seepage/jportiz/nmsba/SlugTide-copy/Gierke_Test/data'  # path of local time data 
data_dir = '..'  # data is in the directory one above this
num_header_rows = 0     # number of header rows to skip  (0 if no header)
nan_string = 'NaT'      # missing value string (e.g., 'NaN')
timezone = 'US/Mountain'  # local timezone
#------------------------------------------------------


# If this script called on its own (with no argument)
#---- (0th arg is script name)
if len(sys.argv)==1:
    print('Enter the file name of the data you want to convert to UTC time:')
    print('   (e.g.,  well_dw4.txt )')
    #  filename = input()            # user enters the filename
    filename = os.path.join(data_dir,input())     # user enters the filename
    print('\nOk. Going to convert this file:')
    print(filename)
# If script called with commandline argument
#  else: filename = sys.argv[1]
else:
    filename = os.path.join(data_dir,sys.argv[1])
    print('\nOk. Going to convert this file:')
    print(filename)
print()
print('  Done.')

if num_header_rows == 0: head = None
else: head = num_header_rows

#-------------------------------------------------- 
#  READ IN THE LOCAL DATA AND MAKE TIMEZONE AWARE
#-------------------------------------------------- 
#  df = pd.read_csv(filename,header=head,na_values=nan_string,delimiter='\t')
# This should automatically deteect the delimiters used in the text file
df = pd.read_csv(filename,header=head,na_values=nan_string,sep=None,engine='python')
df.local_time = pd.to_datetime(df[0]) # Assumes 1st column is time

#---- Make the times tz-aware
df.local_time = df.local_time.dt.tz_localize(timezone, ambiguous='NaT', nonexistent='NaT') #<-- works

#-------------------------------------------------- 
#  CONVERT LOCAL TIME TO UTC
#-------------------------------------------------- 
#  print(df.local_time)
print(df.local_time.dt.tz_convert("UTC"))
df.utc_time = df.local_time.dt.tz_convert("UTC")

new_df = df; new_df.iloc[:,0] = df.utc_time

#-------------------------------------------------- 
#  EXPORT AS .csv 
#-------------------------------------------------- 
outfilename = os.path.basename(filename[:-4])+'_utc.csv'
print('Saving UTC water level data as:')
print('   {}'.format(outfilename))
new_df.to_csv(outfilename, header=False, index=False)
#  new_df.to_csv(outfilename, header=False, index=False, na_rep=nan_string)
print('Done.')
