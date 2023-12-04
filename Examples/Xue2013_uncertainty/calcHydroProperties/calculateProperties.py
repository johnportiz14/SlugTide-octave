'''
Calculate S and T given the observed amplitude response and phase shift for time series data.
Follows the assumptions of Hsieh1987:
    - 2D, isotropic, homogeneous, laterally extensive aquifer
Other relevant literature:
    Xue2013
    Xue2016


Run as:
    python calculateProperties.py -w <wellName>
to run only a specific well listed in ``wellInfo`` file.

Run as:
    python calculateProperties.py -s
to use standard deviation in uncertainty propagation.

Run as:
    python calculateProperties.py -c
to use 95% confidence intervals in uncertainty propagation.

Can also run as:
    python calculateProperties.py -sw <wellName>
    or
    python calculateProperties.py -cw <wellName>
to combine options.

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
import pandas as pd
import argparse
#-----------------------------------------------------
#MATPLOTLIBRC PLOTTING PARAMETERS
# Load up sansmath so that math --> helvetica font
# Also need to tell tex to turn on sansmath package
plt.rcParams['text.latex.preamble'] = r'\usepackage{sansmath}'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.labelweight']=u'normal'
#-----------------------------------------------------

def mc_sim(userParamFileName='UserParam.m'):
    '''
    Check UserParam.m file to see if performing uncertainty quantification
    using MC simulation.
    '''
    filename = os.path.join('..', userParamFileName)
    # Check for MCMC variable and whether it is 0 (False) or 1 (True)
    with open(filename,'r') as f:
        mc_ln = False
        for ln in f:
            if ln.startswith('MC'):
            #  if ln.startswith('MCMC'):
                mc_ln = ln
    mc_val = mc_ln.split('=')[1]
    mc_val = int(mc_val.split(';')[0])
    if mc_val == 0:
        mc_run = False
    elif mc_val == 1:
        mc_run = True
    else: print('MC run not correctly specified in UserParam.m')
    return mc_run

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


def optimize_SandT(pars, amp, phase_shift, period, r_case, r_well, b):
    '''
    ----------------------
    Objective function to be minimized.
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

def invert_SandT(amp, phase_shift, period, r_case, r_well, b, optim_routine='differential_evolution',x0=[]):
    #  SS_array=[]
    #  TT_array=[]
    #  success_array=[]
    if 'diff' in optim_routine:
        print('Beginning differential evolution..')
        # Set guess bounds for log(S) and log(T)
        bnds = ( (-6, -1),       # log bounds for S
                 (-8, -1), )     # log bounds for T
        #  results = differential_evolution(optimize_SandT, bounds = bnds, args=(amp, phase_shift, period, r_case, r_well, b), tol=1e-10, atol=1e-10)
        #  results = differential_evolution(optimize_SandT, bounds = bnds, args=(amp, phase_shift, period, r_case, r_well, b), tol=0.01, workers=8)
        results = differential_evolution(optimize_SandT, bounds = bnds, args=(amp, phase_shift, period, r_case, r_well, b), tol=1e-14, workers=8)
        best_S = 10**results.x[0]
        best_T = 10**results.x[1]
        # TEST
        w = calc_w(period)
        E = calc_E(w,best_S,best_T, r_case, r_well)
        F = calc_F(w,best_S,best_T, r_case, r_well)
        A_test = (E**2 + F**2)**(-0.5)
        A_error = (- A_test)**2
        A_g = calc_A(E,F)                   #calculate resulting ampitude response
        phase_shift_g = np.degrees(calc_eta( E,F )) #calculate resulting phase shift
        #
        print('\n\nFor the parameters:')
        print('  Ampl. resp.: {} '.format(amp))
        print('  Phase shift: {}°'.format(phase_shift))
        print('----')
        #  print('  A_guess    : {} '.format(A_g))
        #  print('  phase_shift_guess  : {}° '.format(phase_shift_g))
        print('--------------------------')
        print('Total Error = {}'.format(results.fun))
        print('after {} evaluations...'.format(results.nfev))
        print('S        = {:.4e}      '.format(best_S))
        print('T        = {:.4e} m^2/s'.format(best_T))
        print('----')
        #  # This only works for the first data point
        #  print('From literature (eyeball estimate):')
        #  print('S_actual ~ {:.4e}      '.format(2.34e-4))
        #  print('T_actual ~ {:.4e} m^2/s'.format(4.30e-6))
        print('Differential evolution successful?')
        print(results.success)
        #  SS_array.append(best_S)
        #  TT_array.append(best_T)
        #  success_array.append(results.success)

    elif 'minim' in optim_routine:
        print('Beginning gradient-based minimization...')
        # Check that initial guess array (x0) is provided
        if len(x0)==0:
            print('Need to provide an initial guess array for [ log(S), log(T) ] ...')
        #---------------------------------------- 
        # MINIMIZE FUNCTION
        #---------------------------------------- 
        #  S_de = results.x[0]
        #  T_de = results.x[1]
        bnds = ( (-6, -1),
                 (-8, -1), )
        #  results = minimize(optimize_SandT, x0=x0, args=(amp, phase_shift, period, r_case, r_well, b), bounds = bnds, tol=0.01) #works well, fast
        results = minimize(optimize_SandT, x0=x0, args=(amp, phase_shift, period, r_case, r_well, b), bounds = bnds, tol=1e-20) #works well, fast
        #  results = minimize(optimize_SandT, x0=x0, args=(amp, phase_shift, period, r_case, r_well, b), bounds = bnds, tol=1e-20) #works well, fast
        best_S = 10**results.x[0]
        best_T = 10**results.x[1]
        print('\nFor the parameters:')
        print('  Ampl. resp.: {} '.format(amp))
        print('  Phase shift: {}°'.format(phase_shift))
        print('--------------------------')
        print('Total Error = {}'.format(results.fun))
        print('after {} evaluations...'.format(results.nfev))
        print('S    = {:.4e}      '.format(best_S))
        print('T    = {:.4e} m^2/s'.format(best_T))
        #  #  t.append(amp_data[i,0])
        #  SS_array.append(best_S)
        #  TT_array.append(best_T)

    else:
        print('Options for ``optim_routine`` arg are "differential_evolution" or "minimize"')
        print('Choose a valid option.')
    #  return SS_array, TT_array, success_array
    return best_S, best_T, results.success

#---------------------------------------- 
#       TESTING
#---------------------------------------- 

if __name__== "__main__":
    #-------------------------------------------------- 
    # Define args for user to choose error to use in Uncertainty Propagation 
    #-------------------------------------------------- 
    descr ='Run the ``calculateProperties.py`` script and (optionally) specify error to use.'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('-s', '--standardDeviation' ,  action='store_true', help='Use standard deviation..')
    parser.add_argument('-c', '--confidenceInterval',  action='store_true', help='Use 95% confidence interval..')
    parser.add_argument('-w', '--well', help='Name of individual Well.')
    args = parser.parse_args()

    # What type of error to use for Uncertainty Propagation?
    if args.standardDeviation:
        err_type = 's'
    elif args.confidenceInterval:
        err_type = 'c'
    else:
        err_type='s'  #default 



    #----User Parameters ------------------------------------------
    tau = 12.421 * 3600. #[s]       period of the M2 constituent
    #r_c = 0.08           #[m] radius of well casing
    #r_w = 0.09           #[m] radius of well
    #b = 400.             #[m] thickness of open interval of well
    rho_water = 1000.    #[kg/m3] density of water
    mu = 1e-3            #[Pa s] fluid dynamic viscosity at 20°C
    # READ INDIVIDUAL WELL PARAMETERS FROM A SEPARATE FILE
    with open('wellInfo') as f:
        d = f.read()
    wellInfo = ast.literal_eval(d)
    wellList = wellInfo.keys()
    #  # IF USER RUNS SCRIPT WITH A SINGLE WELL ARGUMENT, ONLY CALCULATE THAT WELL
    #  #---- Must be a well name that exists in wellInfo file
    #  if len(sys.argv)>1: wellList = [sys.argv[1]]
    # Run individual well only?
    if args.well:
        if args.well in wellList:
            wellList = [args.well]
        else:
            print('Well name must be provided as a key in wellInfo file.')
            print('\nCheck spelling and that name exists.')
            sys.exit()
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
    #    prefix = '../Well_WFSD-1_1'
        ampl_data = np.genfromtxt(prefix+'_amp.csv', delimiter=',',missing_values='NaN')
        eta_data = np.genfromtxt(prefix+'_pha.csv', delimiter=',',missing_values='NaN')

        # Read in CI stats from MC simulations (if on)
        if mc_sim()==True:
            if err_type=='s':
                # For now, read in mean amplitude
                ampl_data_err = np.genfromtxt(prefix+'_amp_stats.csv', delimiter=',',missing_values='NaN')
                eta_data_err  = np.genfromtxt(prefix+'_pha_stats.csv', delimiter=',',missing_values='NaN')
            elif err_type=='c':
                # For now, read in mean amplitude
                #  ampl_data_ci = np.genfromtxt(prefix+'_amp_stats.csv', delimiter=',',missing_values='NaN')
                ampl_data_err = np.genfromtxt(prefix+'_amp_ci95.csv', delimiter=',',missing_values='NaN')
                eta_data_err  = np.genfromtxt(prefix+'_pha_ci95.csv', delimiter=',',missing_values='NaN')


        # Solve for these variables 
        S_array = []
        T_array = []

        t = []  #make new time array because some ampl and eta have nans success_array = []  # record when diff evolution algorithm succeeds


        #---------------------------------------------------------------------------
        # start of time loop 
        #---------------------------------------------------------------------------
        if mc_sim()==False:
            # Calculate S and T for each amplitude-response/phase-shift pair in time
            for i in range(len(ampl_data)):
                # If the value is a NaN, skip it and go to next row
                if np.isnan(ampl_data[i,1]) == True: continue
                # Get the Amplitude and Phase at current time 
                ampl = ampl_data[i,1]      #[1/m]
                # Reciprocal of ampl so units are same as in Xue2016 (m/str)
                ampl = 1/ampl              #[m/str]
                eta = eta_data[i,1]        #[°]
                t.append(ampl_data[i,0])   #[d]

                # Calculate Mean S and T  
                S,T,success = invert_SandT(ampl,eta, tau,r_c,r_w, b, optim_routine='differential_evolution')
                S_array.append(S); T_array.append(T); success_array.append(success)
        #---------------------------------------------------------------------------
        # end of time loop 
        #---------------------------------------------------------------------------

        #---------------------------------------------------------------------------
        # start of CI time loop 
        #---------------------------------------------------------------------------
        # Calculate S and T and 95% Confidence Interval
        if mc_sim()==True:
            i_true = -1
            t_ci = [] #times for CIs can be slightly different
            S_array_mean=[];T_array_mean=[]
            S_array_err_L=[];S_array_err_U=[]
            T_array_err_L=[];T_array_err_U=[]
            success_array_err_L=[];success_array_err_U=[]
            for i in range(len(ampl_data)):
                # If the value is a NaN, skip it and go to next row
                if np.isnan(ampl_data_err[i,1]) == True: continue
                i_true+=1
                print(i)
                #  if np.isnan(eta_data_err[i,1]) == True: continue
                t_ci.append(eta_data_err[i,0])
                if err_type=='s':
                    eta_mean  = eta_data_err[i,1]
                    eta_err_L  = eta_mean-eta_data_err[i,2]
                    eta_err_U  = eta_mean+eta_data_err[i,2]
                elif err_type=='c':
                    eta_mean  = eta_data_err[i,1]
                    eta_err_L  = eta_data_err[i,2]
                    eta_err_U  = eta_data_err[i,3]
                # Reciprocal of ampl so units are same as in Xue2016 (m/str)
                ampl_mean = 1./ampl_data_err[i,1]  #[m/str] use MEAN amplitude for now

                # Calculate S and T mean 
                print('Mean')
                print()
                #  # Use differential evolution
                S, T, success = invert_SandT(ampl_mean,eta_mean, tau,r_c,r_w, b, optim_routine='differential_evolution')
                #  #  # Use gradient-based method 
                #  S, T, success = invert_SandT(ampl_mean,eta_mean, tau,r_c,r_w, b, optim_routine='minimize',x0=(-4.,-6))
                S_array_mean.append(S); T_array_mean.append(T)
#  
                # Calculate S and T at LOWER 95% CI
                print('Lower 95% CI')
                print()
                #  # Use differential evolution
                #  S, T, success = invert_SandT(ampl_mean,eta_err_L, tau,r_c,r_w, b, optim_routine='differential_evolution')
                # Use gradient-method 
                # Use mean S and T as initial guess array (x0)
                #  x0 = ( np.log10(np.mean(S_array)), np.log10(np.mean(T_array)) )
                x0 = ( np.log10(S_array_mean[i_true]), np.log10(T_array_mean[i_true]) )
                S, T, success = invert_SandT(ampl_mean,eta_err_L, tau,r_c,r_w, b, optim_routine='minimize',x0=x0)
                S_array_err_L.append(S); T_array_err_L.append(T); success_array_err_L.append(success)

                # Calculate S and T at LOWER 95% CI
                print('Upper 95% CI')
                print()
                #  # Use differential evolution
                #  S, T, success = invert_SandT(ampl_mean,eta_err_L, tau,r_c,r_w, b, optim_routine='differential_evolution')
                # Use gradient-method 
                S, T, success = invert_SandT(ampl_mean,eta_err_U, tau,r_c,r_w, b, optim_routine='minimize',x0=x0)
                S_array_err_U.append(S); T_array_err_U.append(T); success_array_err_U.append(success)

        #---------------------------------------------------------------------------
        # end of time loop 
        #---------------------------------------------------------------------------

        #  scount = 0
        #  for s in success_array:
            #  if s: scount+=1
        #  print('success % = {:.3f}%'.format(scount/len(success_array)*100))
        #  S_stdev = np.std(S_array)
        #  print('st dev of S   = {}'.format(S_stdev))
        #  print('min(S)  = {:e}'.format(min(S_array)))
        #  print('max(S)  = {:e}'.format(max(S_array)))

        #------------------------------------------------------------ 
        # WRITE S AND T TIME SERIES DATA TO CSV FILE
        #------------------------------------------------------------ 

        #------------------------------------------------------------ 
        #---- Storage Coefficient Output
        #------------------------------------------------------------ 
        if mc_sim()==False:
            colData = np.column_stack( (t,S_array) )
            colNames = ['Time', 'Storativity_[-]']
            filename = os.path.join('output',os.path.basename(prefix)+'_storativity.csv')
        else:
            colData = np.column_stack( (t_ci,S_array_mean, S_array_err_L, S_array_err_U) )
            if err_type=='s':
                colNames = ['Time', 'Storativity_[-]', 'Lower', 'Upper']
            elif err_type=='c':
                colNames = ['Time', 'Storativity_[-]', 'CI_95_Lower', 'CI_95_Upper']
            filename = os.path.join('output',os.path.basename(prefix)+'_storativity_uncertainty.csv')
        # Turn into DataFrame
        df = pd.DataFrame(colData, columns=colNames)
        df.to_csv(filename, sep=',', index=False)
        #------------------------------------------------------------ 
        #---- Transmissivity Output 
        #------------------------------------------------------------ 
        if mc_sim()==False:
            colData = np.column_stack( (t,T_array) )
            colNames = ['Time', 'Transmissivity_[m^2/s]']
            filename = os.path.join('output',os.path.basename(prefix)+'_transmissivity.csv')
        else:
            colData = np.column_stack( (t_ci,T_array_mean, T_array_err_L, T_array_err_U) )
            if err_type=='s':
                colNames = ['Time', 'Transmissivity_[m^2/s]', 'Lower', 'Upper']
            elif err_type=='c':
                colData = np.column_stack( (t_ci,T_array_mean, T_array_err_L, T_array_err_U) )
                colNames = ['Time', 'Transmissivity_[m^2/s]', 'CI_95_Lower', 'CI_95_Upper']
            filename = os.path.join('output',os.path.basename(prefix)+'_transmissivity_uncertainty.csv')
        # Turn into DataFrame
        df = pd.DataFrame(colData, columns=colNames)
        df.to_csv(filename, sep=',', index=False)
        #------------------------------------------------------------ 
        #---- Permeability Output 
        #------------------------------------------------------------ 
        if mc_sim()==False:
            # Calculate permeability (k) from T
            k_array = []
            for T in T_array: k_array.append(calc_perm(T,b=b))
            colData = np.column_stack( (t,k_array) )
            colNames = ['Time', 'Permeability_[m^2]']
            filename = os.path.join('output',os.path.basename(prefix)+'_permeability.csv')
        else:
            # Calculate permeability (k) from T
            k_array = []; k_l=[]; k_u=[]
            for i,T in enumerate(T_array_mean):
                k_array.append(calc_perm(T,b=b))
                k_l.append(calc_perm( T_array_err_L[i], b=b ) )
                k_u.append(calc_perm( T_array_err_U[i], b=b ) )
            colData = np.column_stack( (t_ci,k_array, k_l, k_u) )
            if err_type=='s':
                colNames = [ 'Time', 'Permeability_[m^2]', 'Lower', 'Upper' ]
            elif err_type=='c':
                colNames = [ 'Time', 'Permeability_[m^2]', 'CI_95_Lower', 'CI_95_Upper' ]
            filename = os.path.join('output',os.path.basename(prefix)+'_permeability_uncertainty.csv')
        # Turn into DataFrame
        df = pd.DataFrame(colData, columns=colNames)
        df.to_csv(filename, sep=',', index=False)


        #------------------------------------------------------------ 
        # CALCULATE TIME-AVERAGED S AND T (AND k) VALUES
        #------------------------------------------------------------ 
        if mc_sim()==False:
            S_avg = np.mean(S_array)
            T_avg = np.mean(T_array)
            k_avg = np.mean(k_array)
        else:
            S_avg = np.mean(S_array_mean)
            T_avg = np.mean(T_array_mean)
            k_avg = np.mean(k_array)
        tname = os.path.join('output','avgProperties')
        with open(tname, 'w') as f:
            f.write('---------------------------\n')
            f.write('Time-Averaged Properties:\n')
            f.write('---------------------------\n')
            f.write('S_avg = {:.4e} [-]\n'.format(S_avg))
            f.write('T_avg = {:.4e} [m2/s]\n'.format(T_avg))
            f.write('k_avg = {:.4e} [m2]\n'.format(k_avg))
        print()
        if 'win' in platform:
            os.system('type '+tname)
        else:
            os.system('cat '+tname)

        if mc_sim() == False:
            #------------------------------------------------------------ 
            # PLOT THE S AND T TIMESERIES 
            #------------------------------------------------------------ 
            tdate = [matlabSerial2python(ti) for ti in t]
            fig, (ax1,ax2) = plt.subplots(2, figsize=(10, 6))
            #---- T
            ax1.scatter(tdate, T_array, marker='o',c='k')
            ax1.set_ylabel(r'Transmissivity [m$^2$ s$^{-1}$]')
            #---- k (on 2nd y-axis)
            ax3 = ax1.twinx()
            mn,mx = ax1.get_ylim()  #get ylim of T plot
            # Don't need to re-plot, just convert the units
            ax3.set_ylim(mn*calc_perm(1,b=b), mx*calc_perm(1,b=b))
            ax3.set_ylabel(r'Permeability [m$^2$]', rotation=270,va='bottom')
            #---- S
            ax2.scatter(tdate, S_array, marker='o',c='k')
            ax2.set_ylabel(r'Storage coefficient [—]')
            ax2.set_xlabel('Date')
            plt.tight_layout()
            plt.savefig('output/hydrogeologicPropertiesPlot.pdf')
            plt.savefig('output/hydrogeologicPropertiesPlot.png')
            plt.close('all')

        else:
            #------------------------------------------------------------ 
            # PLOT THE S AND T TIMESERIES 
            #------------------------------------------------------------ 
            tdate = [matlabSerial2python(ti) for ti in t_ci]
            fig, (ax1,ax2) = plt.subplots(2, figsize=(10, 6))
            #---- T
            ax1.grid(which='major',axis='x')
            #  ax1.scatter(tdate, T_array, marker='o',c='k')
            #  ax1.scatter(tdate, T_array_mean, marker='o',c='k')
            #  ax1.fill_between(tdate, T_array_err_L, T_array_err_U, alpha=0.5)
            errL = abs( np.array(T_array_mean) - np.array(T_array_err_L) )
            errU = abs( np.array(T_array_mean) - np.array(T_array_err_U) )
            ax1.errorbar(tdate, T_array_mean, yerr=(errL,errU), ls='',marker='.',markersize=5,color='k')
            #  ax1.errorbar(tdate, T_array_mean, yerr=None, ls='',marker='.',color='k')
            ax1.set_ylabel(r'Transmissivity [m$^2$ s$^{-1}$]')
            #---- k (on 2nd y-axis)
            ax3 = ax1.twinx()
            ax3.grid(which='major',axis='both')
            #  ax3.grid(which='major')
            mn,mx = ax1.get_ylim()  #get ylim of T plot
            # Don't need to re-plot, just convert the units
            ax3.set_ylim(mn*calc_perm(1,b=b), mx*calc_perm(1,b=b))
            ax3.set_ylabel(r'Permeability [m$^2$]', rotation=270,va='bottom')
            #---- S
            errL = abs( np.array(S_array_mean) - np.array(S_array_err_L) )
            errU = abs( np.array(S_array_mean) - np.array(S_array_err_U) )
            ax2.errorbar(tdate, S_array_mean, yerr=(errL,errU), ls='',marker='.',markersize=5,color='k',)
            ax2.set_ylabel(r'Storage coefficient [—]')
            ax2.set_xlabel('Date')
            ax2.grid(which='major')
            #
            plt.tight_layout()
            plt.savefig('output/hydrogeologicPropertiesPlot_uncertainty.pdf')
            plt.savefig('output/hydrogeologicPropertiesPlot_uncertainty.png')
            plt.close('all')




