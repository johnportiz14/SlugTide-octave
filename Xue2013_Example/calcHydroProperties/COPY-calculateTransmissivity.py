'''
Calculate S and T
'''
import os,sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ker,kei,kn,kv
sys.path.append('/project/gas_seepage/jportiz/scripts')
from tools import find_nearest
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import csv

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

def calc_eta(E,F):
    '''(Eqn 16) Phase shift '''
    eta = -np.arctan(F/E)
    return eta

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

def optimize_SandT(pars, amp, phase_shift, period, r_case, r_well ):
    '''
    SLIGHTLY DIFFERENT WAY
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
    '''
    stor = 10**pars[0]  # convert from log(S) to S
    tran = 10**pars[1]  # convert from log(T) to T

    #---- Calculate E_guess and F_guess using the guessed values of S and T
    w = calc_w(period)
    E_g = calc_E(w, stor, tran, r_case, r_well)
    F_g = calc_F(w, stor, tran, r_case, r_well)
    #---- Calculate resulting phase shift and amplitude response using guesses
    eta_g = np.degrees(calc_eta( E_g,F_g )) #calculate resulting phase shift
    A_g = calc_A(E_g,F_g)                   #calculate resulting amplitude response
    #---- Calculate the error between guess and true 
    err_amp = ( amp - A_g )**2         # <-- best so far
    err_pha = (phase_shift - eta_g)**2 # <-- best so far
    #  err_amp = ( amp - A_g )**2 / amp
    #  err_amp = ( amp - A_g )**2 / amp
    #  err_pha = abs((phase_shift - eta_g) / phase_shift)
    #  err_amp = A_g / amp - 1
    #  err_pha = eta_g / phase_shift - 1
    #  err_amp = abs(A_g-amp) / amp - 1
    #  err_pha = abs(eta_g-phase_shift) / phase_shift - 1
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
    #  print('=====')
    #  print('TOTAL ERROR = {}'.format(err_total))
    #  print('    err_amp = {}'.format(err_amp))
    #  print('    err_pha = {}'.format(err_pha))
    #  print('=====')
    #  print('S_guess = {:e}'.format(stor))
    #  print('T_guess = {:e}'.format(tran))
    #  print('-----')
    #  print('A_true = {}'.format(amp))
    #  print('A_guess = {}'.format(A_g))
    #  print()
    #  print('eta_true = {}'.format(phase_shift))
    #  print('eta_guess = {}'.format(eta_g))
    #  print('-------------------------------------')
    return err_total


#---------------------------------------- 
#       TESTING
#---------------------------------------- 

if __name__== "__main__":
    # Average properties for Xue2013 paper
    #  ang_freq_deghr = 28.984    #[°/hr]  angular freq of M2 constituent
    #  ang_freq_degs  = ang_freq_deghr/3600.    #[°/s]  angular freq of M2 constituent
    tau = 12.421 * 3600. #[s]       period of the M2 constituent
    freq_Hz = 1/tau      #[1/s]     freq   of the M2 constituent
    ang_freq = freq_Hz*2.*np.pi #   angular freq  ''     ''
    r_w = 0.08           #[m] radius of well
    r_c = 0.09           #[m] radius of well casing
    user_eta = -20.                  #[°] user-specied phase shift




    print('Beginning differential evolution..')
    # Set fixed vals
    user_eta = -20. #[°]   phase shift
    ampl = 6.34e-7   #[1/m] amplitude response
    #  angular_freq = 1.9324/86400. #[1/s]  ????
    # Set guess bounds for log(S) and log(T)
    bnds = ( (-10, -3.5),    # log bounds for S
             (-8, -4), )     # log bounds for T

    results = differential_evolution(optimize_SandT, bounds = bnds, args=(ampl, user_eta, tau, r_c, r_w), tol=1.e-10)
    print('Differential evolution successful?')
    print(results.success)
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
    print('  Phase shift: {}°'.format(user_eta))
    print('--------------------------')
    print('Total Error = {}'.format(results.fun))
    print('after {} evaluations...'.format(results.nfev))
    print('S    = {:.4e}      '.format(best_S))
    print('T    = {:.4e} m^2/s'.format(best_T))


# USEFUL BIT OF CODE TO WRITE VALS TO CSV FILE #  with open('output/SvsT.csv','w') as f:
        #  writer = csv.writer(f, delimiter=',')
        #  S = []; T = []
        #  for user_S in userS_array:
            #  # These bounds only used in differential evolution
            #  min_dimT = 1*r_c**2/tau    #minimum dimensionless T to look over (1)
            #  max_dimT = 1000*r_c**2/tau #maximum       ''      ''''  ''   ''  (1000)
            #  bnds = ( (np.log10(min_dimT), np.log10(max_dimT)), )
#  
            #  #---- Use either minimize or differential evolution (differential evolution better if guess is bad)
            #  #-- Minimize --
            #  #  results = minimize(optimize_T, x0=(maxT-minT)/2, args=(user_S, user_eta, tau), tol=0.001) #works well, fast
            #  #-- Differential Evolution --
            #  results = differential_evolution(optimize_T, bounds = bnds, args=(user_S, user_eta, tau), tol=0.001) #works well
            #  #---- Get the resulting T that minimizes the difference in phase shift
            #  best_T = results.x[0]           #returns best value chosen for log(T)
            #  print('\nFor the parameters:')
            #  print('  Phase shift: {} deg'.format(user_eta))
            #  print('  Storativity: {:.2e}'.format(user_S))
            #  print('--------------------------')
            #  print('dimT = {:.2f}'.format(10**best_T*tau/r_c**2))
            #  print('T    = {:.4e} m^2/s'.format(10**best_T))
            #  #---- Write S,T to a two-column .csv file
            #  #  S.append(user_S); T.append(best_T)
            #  writer.writerow( (user_S,10**best_T) )
#  
