'''
Analyze the phase response calculations to locate the times of big phase jumps
(may indicate unrecorded timezone/daylight savings time changes).
Requires you to run SlugTide successfully first so that this script can analyze
the phase shift response.

Can be run several ways:
------------------------
>>> python find_data_jumps.py

--> will run the script for all the wells in the
``wellInfo`` file

>>> python find_data_jumps.py <well>

--> will run the script for just the well specified in
the commandline argument. <well> must be a primary
well key in the ``wellInfo`` file.

'''
import os,sys
import pandas as pd
import numpy as np
np.seterr(all='ignore')  #suppress the divide by zero warnings
import matplotlib
import matplotlib.pyplot as plt
from os.path import join
import scipy.fft
from scipy import signal
from scipy.signal import blackman
from scipy import interpolate
import time
from datetime import datetime, timedelta
import glob
import ast

#-----------------------------------------------------
#MATPLOTLIBRC PLOTTING PARAMETERS
# Load up sansmath so that math --> helvetica font
# Also need to tell tex to turn on sansmath package
plt.rcParams['text.latex.preamble'] = r'\usepackage{sansmath}'#,
#    r'\sansmath']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.labelweight']=u'normal'
#-----------------------------------------------------

def floatHourToTime(fh):
    hours, hourSeconds = divmod(fh, 1)
    minutes,seconds = divmod(hourSeconds * 60, 1)
    return (
        int(hours),
        int(minutes),
        int(seconds * 60),
    )

def excelSerial2python(serial_date):
    "Convert an Excel serial date to a Python datetime object"
    # Get Excel serial dates start at  1/1/1900
    dt = datetime.fromordinal(datetime(1900, 1, 1).toordinal() + int(serial_date) - 2)
    tt = dt.timetuple()
    hour, minute, second = floatHourToTime(serial_date % 1)
    dt = dt.replace(hour=hour, minute=minute, second=second)
    return dt

def matlabSerial2python(serial_date):
    "Convert a Matlab serial date to a Python datetime object"
    # Matlab's date0 is 1899-12-30
    date0 = 693960
    # Subtract this from the serial datenum to get Excel serial 
    serial_date_E = serial_date-date0
    dt = excelSerial2python(serial_date_E)
    return dt

if __name__== "__main__":
    #----User Parameters ------------------------------------------
    # READ INDIVIDUAL WELL PARAMETERS FROM A SEPARATE FILE ('wellInfo')
    with open('wellInfo') as f:
        d = f.read()
    wellInfo = ast.literal_eval(d)
    wellList = wellInfo.keys()
    # IF USER RUNS SCRIPT WITH A SINGLE WELL ARGUMENT, ONLY CALCULATE THAT WELL
    #---- Must be a well name that exists as a primary key in  wellInfo file
    print()
    if len(sys.argv)>1:
        wellList = [sys.argv[1]]
        print('Running only for well {}...'.format(wellList[0]))
    else: print('Running for all wells in ``wellInfo`` file...')
    #--------------------------------------------------------------
    
    # Make 'output/' directory if doesn't exist
    if not os.path.exists(os.path.join(os.getcwd(),'output')):
        os.makedirs('output')
    subdir = join('output','spectra')
    if not os.path.exists(os.path.join(os.getcwd(),subdir)):
        os.makedirs(subdir)
    # Detect user operating system
    platform = sys.platform; platform = platform.lower()


    datadir = '..'  #assumes the SlugTide output is all in the directory 1 above

    for well in wellList:
        print()
        print('Well : '+well)
        #  prefix = '../{}'.format(wellInfo[well]['prefix'])
        prefix = join(datadir, wellInfo[well]['prefix'])
        
    #prefix = '../Well-2_1'

    # Read in the phase shift data
    eta_data = np.genfromtxt(prefix+'_pha.csv', delimiter=',',missing_values='NaN')

    #-----------------------------------------------------
    #   FIND BIG JUMPS 
    #-----------------------------------------------------
    # Note, this method prob doesn't work if the delta t is smaller
    #---- should try smoothing/mving average first, then derivative
    diff = np.diff(eta_data[:,1])
    #Big jumps
    thresh = 5. #[°] phase shift jump threshold
    #  big_jumps_ind = np.where(diff > thresh)[0] + 1      #the +1 accounts for the diff shortening(?)
    big_jumps_ind = np.where( abs(diff) > thresh )[0] + 1      #the +1 accounts for the diff shortening(?)

    # Get those times
    big_jumps_serial = eta_data[:,0][big_jumps_ind]
    big_jumps_pytime = [matlabSerial2python(bjs) for bjs in big_jumps_serial]
    tname = os.path.join('output','phase_jumps_list_{}'.format(os.path.basename(prefix)))
    with open(tname, 'w') as f:
        f.write('Detected {} big jumps\n'.format(len(big_jumps_pytime)))
        f.write('  - threshold = {}°\n'.format(thresh))
        f.write('------------------------------------------------\n')
        [f.write(str(bjp)+'\n') for bjp in big_jumps_pytime]
    if ('win' in platform)&('darwin' not in platform):
        os.system('type '+tname)
    else:
        os.system('cat '+tname)
        #print()
#    print('Detected {} big jumps'.format(len(big_jumps_pytime)))
#    print('  - threshold = {}°\n'.format(thresh))
#    [print(bjp) for bjp in big_jumps_pytime]
#    print()
    print('  Done.')

    #-----------------------------------------------------
    #   PLOT 
    #-----------------------------------------------------
    print('Plotting...')
    dates = [matlabSerial2python(s) for s in eta_data[:,0]]
    fig, ax = plt.subplots(1,figsize=(10,4))
    ax.plot(dates, eta_data[:,1],ls = '', marker='.')
    #---- Highlight the big jumps
    win = 1  #window on either side to highlight
    
    for r in big_jumps_pytime:
        #  print(r)
        rp = r+timedelta(days=win)
        rm = r-timedelta(days=win)
        ax.axvspan(rm,rp, color='yellow',alpha=0.75)
    #---- Properties
    ax.hlines(0., ax.get_xlim()[0],ax.get_xlim()[-1],ls='--',color='grey')
    ax.set_ylabel(r'Phase shift [$^\circ$]')
    plt.savefig('output/phase_jumps_{}.pdf'.format(os.path.basename(prefix)))
    plt.close('all')
    print('  Done.')
    print('\nDone with everything.')
    
   



