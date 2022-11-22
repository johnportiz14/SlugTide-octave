'''
Calculate S and T given the observed amplitude response and phase shift for time series data.
Follows the assumptions of Hsieh1987:
    - 2D, isotropic, homogeneous, laterally extensive aquifer
Other relevant literature:
    Xue2013
    Xue2016
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
plt.rcParams['text.latex.preamble'] = r'\usepackage{sansmath}'
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


#---------------------------------------- 
#       TESTING
#---------------------------------------- 

if __name__== "__main__":
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
    #    prefix = '../Well_WFSD-1_1'
        ampl_data = np.genfromtxt(prefix+'_amp.csv', delimiter=',',missing_values='NaN')
        eta_data = np.genfromtxt(prefix+'_pha.csv', delimiter=',',missing_values='NaN') 



        # Solve for these variables 
        S_array = []
        T_array = []

        t = []  #make new time array because some ampl and eta have nans
        success_array = []  # record when diff evolution algorithm succeeds
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
            results = differential_evolution(optimize_SandT, bounds = bnds, args=(ampl, eta, tau, r_c, r_w, b), tol=1e-10, atol=1e-10)
            best_S = 10**results.x[0]
            best_T = 10**results.x[1]
            # TEST
            w = calc_w(tau)
            E = calc_E(w,best_S,best_T, r_c, r_w)
            F = calc_F(w,best_S,best_T, r_c, r_w)
            A_test = (E**2 + F**2)**(-0.5)
            A_error = (ampl - A_test)**2
            A_g = calc_A(E,F)                   #calculate resulting amplitude response
            eta_g = np.degrees(calc_eta( E,F )) #calculate resulting phase shift
            #
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
            #  # This only works for the first data point
            #  print('From literature (eyeball estimate):')
            #  print('S_actual ~ {:.4e}      '.format(2.34e-4))
            #  print('T_actual ~ {:.4e} m^2/s'.format(4.30e-6))


            print('Differential evolution successful?')
            print(results.success)
            success_array.append(results.success)
            #  #---------------------------------------- 
            #  #---------------------------------------- 
            #  # MINIMIZE FUNCTION
            #  #---------------------------------------- 
            #  # Try running optimize after differential evolution to fine-tune the guesses...
            #  print()
            #  print('Now beginning minimization...')
            #  S_de = results.x[0]
            #  T_de = results.x[1]
            #  # new bounds for the minimize function
            #  bnds = ( (-6, -4),
                     #  (-8, -4), )
            #  results = minimize(optimize_SandT, x0=(S_de,T_de), args=(ampl, user_eta, tau, r_c, r_w), bounds = bnds, tol=1e-20) #works well, fast
            best_S = 10**results.x[0]
            best_T = 10**results.x[1]
            print('\nFor the parameters:')
            print('  Ampl. resp.: {} '.format(ampl))
            print('  Phase shift: {}°'.format(eta))
            print('--------------------------')
            print('Total Error = {}'.format(results.fun))
            print('after {} evaluations...'.format(results.nfev))
            print('S    = {:.4e}      '.format(best_S))
            print('T    = {:.4e} m^2/s'.format(best_T))

            t.append(ampl_data[i,0])
            S_array.append(best_S)
            T_array.append(best_T)

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

        #------------------------------------------------------------ 
        # CALCULATE AVERAGE S AND T (AND k) VALUES
        #------------------------------------------------------------ 
        S_avg = np.mean(S_array)
        T_avg = np.mean(T_array)
        k_avg = np.mean(k_array)
        tname = os.path.join('output','avgProperties')
        with open(tname, 'w') as f:
            f.write('------------------------\n')
            f.write('Average Properties:\n')
            f.write('------------------------\n')
            f.write('S_avg = {:.4e} [-]\n'.format(S_avg))
            f.write('T_avg = {:.4e} [m2/s]\n'.format(T_avg))
            f.write('k_avg = {:.4e} [m2]\n'.format(k_avg))
        print()
        if 'win' in platform:
            os.system('type '+tname)
        else:
            os.system('cat '+tname)  

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


