'''
Calculate S and T given the observed amplitude response and phase shift for time series data.
----------------------
Horizontal Flow Model:
----------------------
Follows the assumptions of Hsieh1987:
    - 2D, isotropic, homogeneous, laterally extensive aquifer
----------------------
Vertical Flow Model:
----------------------
Follows the assumptions of Wang2000:
    - vertical flow
Other relevant literature:
    Xue2013
    Xue2016

Can be run several ways:
------------------------
>>> python calculateProperties.py

--> will run the script for all the wells in the
``wellInfo`` file

>>> python calculateProperties.py <well>

--> will run the script for just the well specified in
the commandline argument. <well> must be a primary
well key in the ``wellInfo`` file.

'''
import os,sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ker,kei,kn,kv
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import csv
import glob
import ast
from datetime import datetime

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

def calc_eta(E,F):
    '''(Eqn 16) Phase shift '''
    phase_shift = -np.arctan(F/E)
    return phase_shift

def calc_A(E,F):
    '''(Eqn 15) Amplitude response '''
    A = (E**2 + F**2 )**(-0.5)
    return A

def calc_w(period):
    ''' Calculate angular frequency given period [rad/timeUnit].'''
    return 2*np.pi/period

def calc_tau(angular_freq):
    ''' Calculate period given angular frequency.'''
    return 2*np.pi/angular_freq

def calc_alpha_w(w,S,T,r_well):
    '''
    (Eqn 10) Alpha_omega

    :param w: angular frequency of the fluctuation [1/s]
    :type w: float
    :param S: storage coefficient [1/m]
    :type S: float
    :param T: transmissivity [m^2/s]
    :type T: float
    :param r_well: radius of the screened/open portion of well [m]
    :type r_well: float
    '''
    a_w = ( (w * S)/T )**0.5 * r_well
    return a_w

def calc_E(w, S, T, r_casing, r_well):
    '''
    (Eqn 17) Calculate E

    :param w: angular frequency of the fluctuation [1/s]
    :type w: float
    :param S: storage coefficient [1/m]
    :type S: float
    :param r_casing: radius of the casing
    :type r: float
    :param T: transmissivity [m^2/s]
    :type T: float
    :param r_well: radius of the screened/open portion of well [m]
    :type r_well: float
    '''
    a_w = calc_alpha_w(w,S,T,r_well)
    E = 1 - ( w * r_casing**2. ) / 2 / T * kei(a_w)
    return E

def calc_F(w, S, T, r_casing, r_well):
    '''
    (Eqn 18) Calculate F

    :param w: angular frequency of the fluctuation [1/s]
    :type w: float
    :param S: storage coefficient [1/m]
    :type S: float
    :param r_casing: radius of the casing
    :type r: float
    :param T: transmissivity [m^2/s]
    :type T: float
    :param r_well: radius of the screened/open portion of well [m]
    :type r_well: float
    '''
    a_w = calc_alpha_w(w,S,T,r_well)
    F = ( w * r_casing**2. ) / 2 / T * ker(a_w)
    return F

def dim_omega(w, r_casing, T):
    return w * r_casing**2 / T

def dim_stor(S, r_casing, r_well):
    return S * r_well**2 / r_casing**2

def calc_perm(T, b,  mu=1e-3, rho=1000, g=9.8):
    '''
    Calculate permeability (k; [m^2]) given the tranmissivity (T; [m^2/s]).

    :param T: transmissivity [m^2/s].
    :type T: float
    :param b: thickness of open interval of the well [m].
    :type b: float
    :param mu: dynamic viscosity of the fluid [Pa s].
    :type mu: float
    :param rho: density of the liquid [kg/m^3].
    :type rho: float
    :param g: gravitational acceleration [m/s^2].
    :type g: float
    '''
    k = mu / (rho*g*b) * T
    return k


def optimize_SandT_horz(pars, amp, phase_shift, period, r_case, r_well, b):
    '''
    ----------------------
    Objective function to be minimized.
    - Horizontal flow model (Hsieh1987)
    - Modify values of S and T simultaneously and compare resulting values
    of amplitude response and phase shift.

    :param pars: input values to be guessed (provided by
    differential_evolution function).
        pars[0] will be Storage coefficient [-]
        pars[1] will be Transmissovity [m2/s]
    :type pars: scipy object (?)
    :param amp: observed amplitude response [1/m].
    :type amp: float
    :param phase_shift: observed phase shift [°].
    :type phase_shift: float
    :param period: period of the tidal forcing [s] (e.g., period of M2 = 12.421 hr --> 44715.6 s).
    :type period: float
    :param r_case: inner radius of well casing [m].
    :type r_case: float
    :param r_well: radius of well [m].
    :type r_well: float
    :param b: well open interval thickness [m].
    :type b: float
    '''
    stor = 10**pars[0]  # convert from log(S) to S
    tran = 10**pars[1]  # convert from log(T) to T
    #  stor = pars[0]  # convert from log(S) to S
    #  tran = pars[1]  # convert from log(T) to T

    # !!!! NEW !!!!
    Ss = stor/b         #[1/m] calculate Specific storage 
    print('S  = {}'.format(stor))
    print('Ss = {}\n'.format(Ss))

    #---- Calculate E_guess and F_guess using the guessed values of S and T
    w = calc_w(period)
    E_g = calc_E(w, stor, tran, r_case, r_well)
    F_g = calc_F(w, stor, tran, r_case, r_well)
    #---- Calculate resulting phase shift and amplitude response using guesses
    eta_g = np.degrees(calc_eta( E_g,F_g )) #calculate resulting phase shift
    print('eta_g  = {}'.format(eta_g))                             #DEBUG
    print('eta_actual = {}'.format(phase_shift))
    A_g = 1/Ss * calc_A(E_g,F_g)                   #calculate resulting amplitude response

    print('A_g   = {:e}'.format(A_g))                       #DEBUG
    #---- Calculate the error between guess and true 
    #  err_amp = ( amp - A_g )**2         # <-- best so far [JPO] TEST
    #  err_amp = ( np.log10(amp) - np.log10(A_g) )**2         #  !!! TEST !!!
    #  print('A_actual = {}'.format(1/Ss*amp))
    print('A_actual = {:e}'.format(amp))
    #  err_amp = ( amp - A_g )**2         # <-- best so far [JPO] TEST    try Xue2016 way 
    err_amp = ( np.log10(amp) - np.log10(A_g) )**2         #  !!! TEST !!!  try Xue2016 way
    #  err_amp = ( 1/Ss*amp - A_g )**2         # <-- best so far [JPO] TEST    try Xue2016 way 
    #  err_amp = ( np.log10(1/Ss*amp) - np.log10(A_g) )**2         #  !!! TEST !!!  try Xue2016 way
    err_pha = (phase_shift - eta_g)**2 # <-- best so far
    err_total = err_amp + err_pha
    #  err_total = (err_amp + err_pha)/2
    #  err_total =  err_pha
    #  #---- Warning if the phase angle is way off
    #  phase_diff = eta_g - phase_shift
    #  phase_tol = 2. #[°] arbitrary max diff in phase shift
    #  if abs(phase_diff)>phase_tol:
        #  print('***** WARNING ************************************')
        #  print('eta_g = {}'.format(eta_g))
        #  print('Calculated phase shift is off by more than {:.2}°.'.format(phase_tol))
        #  print('**************************************************')
    print('=====')
    print('TOTAL ERROR = {}'.format(err_total))
    print('    err_amp = {}'.format(err_amp))
    print('    err_pha = {}'.format(err_pha))
    print('=====')
    return err_total

def optimize_SandT_vert(pars, amp, phase_shift, period, r_case, r_well, b, z):
    '''
    ----------------------
    Objective function to be minimized.
    - Vertical flow model (Wang200)
        - eqn 6,7 in Xue2016
    - Modify values of S and T simultaneously and compare resulting values
    of amplitude response and phase shift.

    :param pars: input values to be guessed (provided by
    differential_evolution function).
        pars[0] will be Storage coefficient [-]
        pars[1] will be Transmissovity [m2/s]
    :type pars: scipy object (?)
    :param amp: observed amplitude response [1/m].
    :type amp: float
    :param phase_shift: observed phase shift [°].
    :type phase_shift: float
    :param period: period of the tidal forcing [s] (e.g., period of M2 = 12.421 hr --> 44715.6 s).
    :type period: float
    :param r_case: inner radius of well casing [m].
    :type r_case: float
    :param r_well: radius of well [m].
    :type r_well: float
    :param b: well open interval thickness [m].
    :type b: float
    :param z: depth from the surface [m].
    :type z: float
    '''
    stor = 10**pars[0]  # convert from log(S) to S
    tran = 10**pars[1]  # convert from log(T) to T

    w = calc_w(period)
    Ss = stor/b         #[1/m] calculate Specific storage 
    D = tran/stor  #[m2/s] hydraulic diffusivity
    kron_delta = np.sqrt((2*D)/w)
    #---- Calculate A_guess and eta_guess using the guessed values of S and T
    A_g = 1/Ss * (1 - 2.*np.exp(-z/kron_delta)*np.cos(z/kron_delta) + np.exp((-2.*z)/kron_delta))**0.5
    eta_g = np.degrees( np.arctan( (np.exp(-z/kron_delta)*np.sin(z/kron_delta)) / (1-np.exp(-z/kron_delta)*np.cos(z/kron_delta)) ) )
    print('eta_g  = {}'.format(eta_g))                             #DEBUG
    print('eta_actual = {}'.format(phase_shift))


    print('A_g   = {:e}'.format(A_g))                       #DEBUG
    #---- Calculate the error between guess and true 
    print('A_actual = {:e}'.format(amp))
    err_amp = ( np.log10(amp) - np.log10(A_g) )**2         #  !!! TEST !!!
    err_pha = (phase_shift - eta_g)**2 # <-- best so far
    err_total = err_amp + err_pha
    print('=====')
    print('TOTAL ERROR = {}'.format(err_total))
    print('    err_amp = {}'.format(err_amp))
    print('    err_pha = {}'.format(err_pha))
    print('=====')
    return err_total



#---------------------------------------- 
#       TESTING
#---------------------------------------- 

if __name__== "__main__":
    #----User Parameters ------------------------------------------
    tau = 12.421 * 3600. #[s]       period of the M2 constituent
    rho_water = 1000.    #[kg/m3] density of water
    mu = 1e-3            #[Pa s] fluid dynamic viscosity at 20°C
    g = 9.81             #[m2/s] gravitational acceleration
    # READ INDIVIDUAL WELL PARAMETERS FROM A SEPARATE FILE
    with open('wellInfo') as f:
        d = f.read()
    wellInfo = ast.literal_eval(d)
    wellList = wellInfo.keys()
    # IF USER RUNS SCRIPT WITH A SINGLE WELL ARGUMENT, ONLY CALCULATE THAT WELL
    #---- Must be a well name that exists in wellInfo file
    if len(sys.argv)>1: wellList = [sys.argv[1]]
    #--------------------------------------------------------------

    # Make 'output/' directory if doesn't exist
    if not os.path.exists(os.path.join(os.getcwd(),'output')):
        os.makedirs('output')

		# Detect user operating system
    platform = sys.platform; platform = platform.lower()

    #--------------------------------------------------------------
    #  CALCULATE HYDROGEOLOGIC PROPERTIES FOR EACH WELL
    #--------------------------------------------------------------
    for well in wellList:
        print('well : '+well)
        prefix = '../{}'.format(wellInfo[well]['prefix'])
        r_c = wellInfo[well]['r_c']
        r_w = wellInfo[well]['r_w']
        b   = wellInfo[well]['b']
        flowModel = wellInfo[well]['flowModel'] #horizontal or vertical flow model?
        if 'vert' in flowModel:
            try: z = wellInfo[well]['z']  #probably only need depth for vert flow
            except KeyError:
                print('\nERROR: you specified a vertical flow model but did not provide a value for "z".')
                break


        # Read the amplitude and phase responses from csv files (skip nans)
        #  prefix = '../Well_WFSD-1_1'
        ampl_data = np.genfromtxt(prefix+'_amp.csv', delimiter=',',missing_values='NaN')
        eta_data = np.genfromtxt(prefix+'_pha.csv', delimiter=',',missing_values='NaN')
        # Get the time bounds for original data (for plotting)
        tmin = ampl_data[0,0]
        tmax = ampl_data[-1,0]


        # Solve for these variables 
        S_array = []
        T_array = []

        t = []  #make new time array because some ampl and eta have nans
        success_array = []  # record when diff evolution algorithm succeeds
        err_array     = []  # record error from optimization routine
        # Calculate S and T for each amplitude-response/phase-shift pair in time
        for i in range(len(ampl_data)):
            # If the value is a NaN, skip it and go to next row
            if np.isnan(ampl_data[i,1]) == True: continue

            ampl = ampl_data[i,1]
            # Invert ampl so units are same as in Xue2016 (m/str)
            ampl = 1/ampl
            eta = eta_data[i,1]

            print('Beginning differential evolution..')
            # Set guess bounds for log(S) and log(T)
            bnds = ( (-6, -1),       # log bounds for S
                     (-8, -1), )     # log bounds for T
            #---------------------------------------------
            # HORIZONTAL OR VERTICAL FLOW MODEL
            #---------------------------------------------
            if 'hor' in flowModel:
                print('Assuming a horizontal flow model...')
                # If eta is positive, don't try to calculate S and T
                #---- Horizontal model range of eta is -90° to 0° 
                #  if eta >0. : continue
                #  if eta >-1. : continue
                if (eta>0. or eta<-90.): continue
                else:
                    results = differential_evolution(optimize_SandT_horz, bounds = bnds, args=(ampl, eta, tau, r_c, r_w, b), tol=1e-10, atol=1e-10)
            elif 'vert' in flowModel:
                print('Assuming a vertical flow model...')
                # If eta is negative, don't try to calculate S and T
                #---- Vertical model range of eta is -45° to -1° 
                #  if eta <0. : continue
                #  if eta <0. : continue
                if (eta<-1. or eta>45.): continue
                else:
                    results = differential_evolution(optimize_SandT_vert, bounds = bnds, args=(ampl, eta, tau, r_c, r_w, b,z), tol=1e-10, atol=1e-10)
                    #  results = differential_evolution(optimize_SandT_vert, bounds = bnds, args=(ampl, eta, tau, r_c, r_w, b,z), tol=1e-20, atol=1e-20,popsize=500)
            else: sys.exit('Did not detect a valid flowModel (default="horizontal".')

            best_S = 10**results.x[0]
            best_T = 10**results.x[1]
            print('\n\nFor the parameters:')
            print('  Ampl. resp.: {} '.format(ampl))
            print('  Phase shift: {}°'.format(eta))
            print('----')
            #  print('  A_guess    : {} '.format(A_g))
            #  print('  eta_guess  : {}° '.format(eta_g))
            print('--------------------------')
            print('Total Error = {}'.format(results.fun))
            print('after {} evaluations...'.format(results.nfev))
            print('S        = {:.4e}      '.format(best_S))
            print('T        = {:.4e} m^2/s'.format(best_T))
            print('----')
            print('Differential evolution successful?')
            print(results.success)
            success_array.append(results.success)

            t.append(ampl_data[i,0])
            S_array.append(best_S)
            T_array.append(best_T)
            err_array.append(results.fun)

        scount = 0
        for s in success_array:
            if s: scount+=1
        print('success % = {:.3f}%'.format(scount/len(success_array)*100))
        S_stdev = np.std(S_array)
        print('st dev of S   = {}'.format(S_stdev))
        print('min(S)  = {:e}'.format(min(S_array)))
        print('max(S)  = {:e}'.format(max(S_array)))

        #------------------------------------------------------------ 
        # WRITE S AND T TIME SERIES DATA TO CSV FILE
        #------------------------------------------------------------ 
        #---- Storage Coefficient Output
        twocols = np.column_stack( (t,S_array) )
        with open(os.path.join('output',os.path.basename(prefix)+'_storativity.csv'), 'w', newline='') as f:
            writer = csv.writer(f,delimiter=',')
            writer.writerows(twocols)
        #---- Transmissivity Output 
        twocols = np.column_stack( (t,T_array) )
        with open(os.path.join('output',os.path.basename(prefix)+'_transmissivity.csv'), 'w', newline='') as f:
            writer = csv.writer(f,delimiter=',')
            writer.writerows(twocols)
        #---- Permeability Output 
        # Calculate permeability (k) from T
        k_array = []
        for T in T_array: k_array.append(calc_perm(T,b=b))
        twocols = np.column_stack( (t,k_array) )
        with open(os.path.join('output',os.path.basename(prefix)+'_permeability.csv'), 'w', newline='') as f:
            writer = csv.writer(f,delimiter=',')
            writer.writerows(twocols)
        #---- Optimization Errors (total error) 
        twocols = np.column_stack( (t,err_array) )
        with open(os.path.join('output',os.path.basename(prefix)+'_error.csv'), 'w', newline='') as f:
            writer = csv.writer(f,delimiter=',')
            writer.writerows(twocols)

        #------------------------------------------------------------ 
        # CALCULATE AVERAGE S AND T (AND k) VALUES
        #------------------------------------------------------------ 
        S_avg = np.mean(S_array)
        S_s_avg = S_avg/b
        T_avg = np.mean(T_array)
        k_avg = np.mean(k_array)
        K_avg = (np.mean(k_array)*rho_water*g)/mu
        tname = os.path.join('output','avgProperties_{}'.format(os.path.basename(prefix)))
        with open(tname, 'w') as f:
            f.write('------------------------\n')
            f.write('Average Properties ({}):\n'.format(os.path.basename(prefix)))
            f.write('------------------------\n')
            f.write('S_avg   = {:.4e} [-]\n'.format(S_avg))
            f.write('S_s_avg = {:.4e} [1/m]\n'.format(S_s_avg))
            f.write('T_avg   = {:.4e} [m2/s]\n'.format(T_avg))
            f.write('k_avg   = {:.4e} [m2]\n'.format(k_avg))
            f.write('K_avg   = {:.4f} [m/d]\n'.format(K_avg*86400.))
            #  f.write('D_avg   = {:.4e} [m2/s]\n'.format(T_avg/S_avg))
            #  f.write('R_avg   = {:.2f} [m]\n'.format(np.sqrt(T_avg/S_avg*tau)))
        print()
        if 'win' in platform:
            os.system('type '+tname)
        else:
            os.system('cat '+tname)      	

        #------------------------------------------------------------ 
        # PLOT THE S AND T TIMESERIES 
        #------------------------------------------------------------ 
        # Convert the serial dates to a datetime format for plotting
        tdate = [matlabSerial2python(ti) for ti in t]
        fig, (ax1,ax2) = plt.subplots(2, figsize=(10, 6), sharex=True)
        #---- T
        #  ax1.scatter(t, T_array, marker='o',c='k')
        ax1.scatter(tdate, T_array, marker='o',c='k')
        ax1.set_ylabel(r'Transmissivity [m$^2$ s$^{-1}$]')
        ax1.set_yscale('log')
        #  ax1.set_xlim(matlabSerial2python(tmin), matlabSerial2python(tmax))
        #---- k (on 2nd y-axis)
        ax3 = ax1.twinx()
        mn,mx = ax1.get_ylim()  #get ylim of T plot
        # Don't need to re-plot, just convert the units
        ax3.set_ylim(mn*calc_perm(1,b=b), mx*calc_perm(1,b=b))
        ax3.set_ylabel(r'Permeability [m$^2$]', rotation=270,va='bottom')
        #---- S
        #  ax2.scatter(t, S_array, marker='o',c='k')
        ax2.scatter(tdate, S_array, marker='o',c='k')
        #---- Properties
        ax2.set_yscale('log')
        ax3.set_yscale('log')
        ax2.set_ylabel(r'Storage coefficient [—]')
        ax2.set_xlabel('Date')
        ax1.set_title('Well {}'.format(well))
        plt.tight_layout()
        plt.savefig('output/hydrogeologicPropertiesPlot_{}.pdf'.format(os.path.basename(prefix)))
        plt.savefig('output/hydrogeologicPropertiesPlot_{}.png'.format(os.path.basename(prefix)))
        plt.close('all')


    # Skip the 'all' plot if only a single well run using commandline argument
    if len(sys.argv)>1:
        #---------------------------------------------
        # PLOT T AND S FOR ALL WELLS ON SINGLE PLOT 
        #---------------------------------------------
        counter = 0
        fig, (ax1,ax2) = plt.subplots(2, figsize=(10, 6))
        for well in wellList:
            counter+=1
            # Read in data that was saved to csv
            prefix = wellInfo[well]['prefix']
            b = wellInfo[well]['b']
            td = np.genfromtxt(os.path.join('output',prefix+'_transmissivity.csv'),delimiter=',')
            sd = np.genfromtxt(os.path.join('output',prefix+'_storativity.csv'),delimiter=',')
            time = td[:,0]
            tdate = [matlabSerial2python(ti) for ti in time]
            T_array = td[:,1]
            S_array = sd[:,1]
            S_s_array = [s/b for s in S_array]
            k_array = [calc_perm(tran,b) for tran in T_array]
            #---- T
            ax1.scatter(tdate, T_array, marker='o',label=well)
            ax1.set_yscale('log')
            #---- k (on 2nd y-axis)
            if counter == len(wellList):
                ax3 = ax1.twinx()
                mn,mx = ax1.get_ylim()  #get ylim of T plot
                # Don't need to re-plot, just convert the units of the axis limits
                ax3.set_ylim(mn*calc_perm(1,b=b), mx*calc_perm(1,b=b))
                #  ax1.scatter(tdate, k_array, marker='o',label=well)
            #---- S
            ax2.scatter(tdate, S_array, marker='o')
            #  ax2.scatter(tdate, S_s_array, marker='o')
        #Plot Properties
        ax1.legend()
        ax1.set_title('All Wells')
        ax2.set_yscale('log')
        ax3.set_yscale('log')
        ax1.set_ylabel(r'Transmissivity [m$^2$/s]')
        ax3.set_ylabel(r'Permeability [m$^2$]', rotation=270,va='bottom')
        #  ax1.set_ylabel(r'Permeability [m$^2$]')
        ax2.set_ylabel(r'Storativity [--]')
        #  ax2.set_ylabel(r'Specific Storage [m$^{-1}$]')
        ax2.set_xlabel('Date')
        #  ax1.set_xlim(right=datetime(2015,12,1,0,0,0)) #add extra space for legend  on right
        ax2.set_xlim(ax1.get_xlim())
        plt.tight_layout()
        plt.savefig('output/hydrogeologicPropertiesPlot_all.pdf')
        plt.savefig('output/hydrogeologicPropertiesPlot_all.png')
        plt.close('all')


