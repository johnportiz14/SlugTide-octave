
% Parameters used in later scripts to calculate and plot the tidal response in a set of wells.
% Inputs below are listed by the script where they are used.

%MasterWell.m
NWells=1;   % Total number of wells

%MCMC ? Uncertainty Propagation
MC=1;     % 0 (False) or 1 (True)
wl_precision=0.006; %[m] Well precision 
N_mc = 1e4; %1e2;%50;  % number of Monte-Carlo simulations to perform


%LoadTides.m
Tide_dt=2/24/60; %[days] tidal sample spacing in units of days for the synthetic tidal model (check SPOTL input)

%BuildWell.m
dt=10/24/60; %[days] time sample spacing of the well data in units of days
tidedate='2010-01-01 00:00:00'; %date of start of synthetic tide 'YYYY-MM-DD hh:mm:ss' format
tcUTC=0; % correction from time of sy nthetic tide to UTC (if processed in local time).
            % may require some preprocessing if Daylight Savings Time is not already corrected for.
tide_dt=dt; % interpolate the tidal model to have the same dt as the data

%BuildWell.m, PlotOriginalWL.m, PlotWellResponse.m
cUTC=-8/24;%7/24; %[days] correction from time of collected data to UTC (if collected in local time). 
            % may require some preprocessing if Daylight Savings Time is not already corrected for.
            % SIGN CONVENTION: if time zone is reported as [UTC+8], use the
            % opposite sign here (time correction TO UTC).
            
%PlotOriginalWL.m
PlotTime1='2010-01-01'; % Date of beginning of plot in 'YYYY-MM-DD' format
% PlotTime2='2010-02-21'; % Date of ending of plot in 'YYYY-MM-DD' format
% PlotTime1='2010-01-01'; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime2='2011-08-06'; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of the dataset timeframe.

% %PlotOriginalWL.m, PlotWellResponse.m
% EQtime='2010-02-26 00:00:00'; % Date of earthquake in 'YYYY-MM-DD hh:mm:ss' format
% %%%%[JPO] attempt to implement multiple earthquakes
EQtime=[
        '2010-02-26 00:00:00';
        '2010-04-06 00:00:00';
        '2010-04-13 00:00:00';
        '2010-12-21 00:00:00';
        '2011-03-11 00:00:00';
        '2011-03-24 00:00:00'
        ]; % Date(s) of earthquake(s) in 'YYYY-MM-DD hh:mm:ss' format 
            % (if there is an earthquake occurence in the data).

%WaterResp_general_interp.m
freq_low=30/24; % 30-hour cycle --> (0.8 cyc/day)
freq_high=10/24; % 10-hour cycle--> (2.4 cyc/day)
t_win=29.6;%29.5307; % days   %%[JPO] Xue2013 used 29.6-day window (this change didn't make huge diff on its own)
n_shift=0.2; %0.1; % percent of each shift of data segment %%[JPO] Xue2013 had 80% overlap b/w segments
per_gap=0.2;%0.1; % percent gap of the time window  %%[JPO] Xue2013 said "﻿If the gap comprises more than 20% of the segment, no further analysis is attempted for the segment"

%PlotWellResponse.m
PlotTime3=PlotTime1; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime4=PlotTime2; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of the dataset timeframe.
    
%avg_amp_phase.m (JPO: COMMENTED OUT FOR NOW)
PlotTime5='2010-01-01'; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime6='2011-08-06'; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of
    % the timeframe during which the user would like to view the average response.