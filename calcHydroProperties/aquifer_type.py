'''
Plot the power spectrum response of filtered, detrended well water
levels so as to perform Rahi analysis to determine aquifer type.

Can be run several ways:
------------------------
>>> python aquifer_type.py

--> will run the script for all the wells in the
``wellInfo`` file

>>> python aquifer_type.py <well>

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
from datetime import datetime
import glob
import ast

#-----------------------------------------------------
#MATPLOTLIBRC PLOTTING PARAMETERS
# Load up sansmath so that math --> helvetica font
# Also need to tell tex to turn on sansmath package
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{sansmath}',
    r'\sansmath']
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


    datadir = '..'  #assumes the SlugTide output is all in the directory 1 above

    for well in wellList:
        print()
        print('Well : '+well)
        #  prefix = '../{}'.format(wellInfo[well]['prefix'])
        prefix = join(datadir, wellInfo[well]['prefix'])
        #---- Read in the filtered, detrended data output by SlugTide
        try:
            wl = pd.read_csv(prefix+'_filtered_detrended.csv',names=['TIME','WL'])
        except FileNotFoundError:
            print('No filtered, detrended data found for well {}...'.format(well))
            print('  skipping...')
            continue
        # Convert time to hours
        wl['TIME'] = wl['TIME']*24. #SlugTide outputs time in days
        dt_hr = np.mean(np.diff(wl['TIME'])) #[hr] sample interval
        #  baro['TIME'] = baro['TIME']*24.

        # dictionary of 5 Major Harmonic Components (with period in hours)
        tidal_dict = {'O1': [25.8193, 'principal lunar'],
                      'K1': [23.9344, 'lunar-solar'],
                      'M2': [12.4206, 'principal lunar'],
                      'S2': [12.0000, 'principal solar'],
                      #  'N2': [12.6583, 'lunar eliptic'],
                     }

        #--------------------------------------------------------------
        # PLOT THE POWER SPECTRUM  
        #--------------------------------------------------------------
        print('Plotting power spectrum...')
        #  http://faculty.jsd.claremont.edu/jmilton/Math_Lab_tool/Labs/Lab9.pdf
        x = wl['WL']
        sample_rate = 1/dt_hr #[sample per hr]  
        #  NFFT = 2**12#4096#256  #may need to tweak this to resolve certain features  <--- RAHI
        NFFT = 2**14#4096#256  #may need to tweak this to resolve certain features  <--- MARIO
        #  Fs = sample_rate/2
        Fs = sample_rate  #for some reason the sample_rate/2 doubles the period...
        power, freqs = matplotlib.mlab.psd(x, NFFT, Fs)
        #---- Save these to an output subdirectory
        print('Saving power and freqs to {}...'.format(subdir))
        twocols = np.column_stack( (freqs, power) )
        np.savetxt(join(subdir,os.path.basename(prefix)+'_spectra.csv'), twocols, delimiter=',')

        #  power, freqs = matplotlib.mlab.psd(x, NFFT, Fs, noverlap=75)
        fig, ax = plt.subplots(figsize=(10,4))
        ax.plot(1/freqs,power)
        #  ax.plot((2*np.pi)/freqs,power)
        ax.set_xlim(0,40)
        ax.set_xlabel('Period [hrs/cycle]')
        ax.set_ylabel('Power')
        # add in dashed lines for major tidal components
        for comp in tidal_dict:
            ax.axvline(x=tidal_dict[comp][0],ls='--',color='k',lw=0.5)
        # Plot names of the tidal components
        for comp in tidal_dict:
            xshift = 1.5#0.17 # how much to shift the label from the line position
            ypos = 0.95*max(power)#7.0    # default ypos
            # Shift to the left...
            if comp in ['O1', 'M2', 'N2']: xpos = tidal_dict[comp][0]+xshift
            # Shift to the right...
            elif comp in ['K1', 'S2']: xpos = tidal_dict[comp][0]-xshift
            if comp in ['N2']: ypos = 0.95*max(power)#2.5 # make N2 label lower than the others
            ax.text(xpos,ypos, comp,ha='center')
        if 'well' not in well.lower(): welltitle = 'Well {}'.format(well)
        else: welltitle = well
        ax.set_title('{}'.format(welltitle))
        fig.tight_layout()
        plt.savefig(join('output','{}_power_spectrum_plot.pdf'.format(os.path.basename(prefix))))
        plt.close('all')
        #  print('Saving plot in {}...'.format(join(os.getcwd(),'output')))
        print('Saving plot in {} directory...'.format(join('output')))
        print('  Done.')

    #--------------------------------------------------------------
    # PLOT ALL WELLS TOGETHER  (IF APPLICABLE)
    #--------------------------------------------------------------
    # Skip the 'all' plot if only a single well run using commandline argument

    if len(sys.argv)==1:
        numplots = len(wellList)
        #  fig, (ax1,ax2) = plt.subplots(numplots) #, figsize=(10, 6))
        #  fig, axs = plt.subplots(numplots) #, figsize=(10, 6))
        fig, axs = plt.subplots(numplots, figsize=(10, 2*numplots))
        fig.subplots_adjust(hspace = .5, wspace=.5)
        counter = 0
        for well in wellList:
            # Read in data that was saved to csv
            prefix = wellInfo[well]['prefix']
            data = np.genfromtxt(join(subdir,os.path.basename(prefix)+'_spectra.csv'),delimiter=',')
            freqs = data[:,0]; power = data[:,1]
            #---- PLOT --------------------------
            axs[counter].plot(1/freqs,power)
            # Are these ylims more extreme?
            axs[counter].set_xlim(0,40)

            axs[counter].set_ylabel('Power')
            # add in dashed lines for major tidal components
            for comp in tidal_dict:
                axs[counter].axvline(x=tidal_dict[comp][0],ls='--',color='k',lw=0.5)
            # Plot names of the tidal components
            for comp in tidal_dict:
                xshift = 1.5#0.17 # how much to shift the label from the line position
                ypos = 0.95*max(power)#7.0    # default ypos
                # Shift to the left...
                if comp in ['O1', 'M2', 'N2']: xpos = tidal_dict[comp][0]+xshift
                # Shift to the right...
                elif comp in ['K1', 'S2']: xpos = tidal_dict[comp][0]-xshift
                if comp in ['N2']: ypos = 0.95*max(power)#2.5 # make N2 label lower than the others
                axs[counter].text(xpos,ypos, comp,ha='center')
            if 'well' not in well.lower(): welltitle = 'Well {}'.format(well)
            else: welltitle = well
            axs[counter].text(0.97,0.85, welltitle, ha='right',transform=axs[counter].transAxes, weight='bold')
            if counter==numplots: axs[counter].set_xlabel('Period [hrs/cycle]')
            counter+=1
        #  axs.set_title('{}'.format(welltitle))
        #---- Make all yLims the same
        miny=1e9; maxy=-1e9  #get most extreme ylims so all subplots can have same limits
        for i in range(numplots):
            ylims = axs[i].get_ylim()
            if ylims[0]<miny: miny = ylims[0]
            if ylims[1]>maxy: maxy = ylims[1]
        for i in range(numplots):
            axs[i].set_ylim(miny,maxy)

        fig.tight_layout()
        plt.savefig(join('output','allWells_power_spectrum_plot.pdf'.format(os.path.basename(prefix))))
        plt.savefig(join('output','allWells_power_spectrum_plot.png'.format(os.path.basename(prefix))))
        plt.close('all')
        #  print('Saving plot in {}...'.format(join(os.getcwd(),'output')))
        print('Saving plot in {} directory...'.format(join('output')))

    #
    print('\nDone with everything.')





