
% Parameters used in later scripts to calculate and plot the tidal response in a set of wells.
% Inputs below are listed by the script where they are used.

%MasterWell.m
NWells=3;   % Total number of wells

%LoadTides.m
Tide_dt=5/24/60; %[days] tidal sample spacing in units of days for the synthetic tidal model (check SPOTL input)

%BuildWell.m
dt=15/24/60; %[days] time sample spacing of the well data in units of days
tidedate='2012-10-01 00:00:00'; %date of start of synthetic tide 'YYYY-MM-DD hh:mm:ss' format
tcUTC=0; % correction from time of synthetic tide to UTC (if processed in local time).
            % may require som preprocessing if Daylight Savings Time is not already corrected for.
            % [default is 0]
tide_dt=dt; % interpolate the tidal model to have the same dt as the data

%BuildWell.m, PlotOriginalWL.m, PlotWellResponse.m
cUTC=7/24; %[days] correction from time of collected data to UTC (if collected in local time). 
            % may require some preprocessing if Daylight Savings Time is not already corrected for.
            % SIGN CONVENTION: if time zone is reported as [UTC+8], use the
            % negative sign here (time correction TO UTC), and vice versa.
            
%PlotOriginalWL.m
PlotTime1='2012-10-01'; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime2='2013-06-30'; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of the dataset timeframe.

% %PlotOriginalWL.m, PlotWellResponse.m
% % EQtime='2003-12-22 12:16:6'; % Date of earthquake in 'YYYY-MM-DD hh:mm:ss' format 
%             % (if there is an earthquake occurence in the data).
% %%%%[JPO] attempt to implement multiple earthquakes capability
% EQtime=[
%         '2010-02-26 00:00:00';
%         '2010-04-06 00:00:00';
%         '2010-04-13 00:00:00';
%         '2010-12-21 00:00:00';
%         '2011-03-11 00:00:00';
%         '2011-03-24 00:00:00'
%         ]; % Date(s) of earthquake(s) in 'YYYY-MM-DD hh:mm:ss' format 
%             % (if there is an earthquake occurence in the data).

%WaterResp_general_interp.m
freq_low=30/24;  %[d] 30-hour cycle --> (0.8 cyc/day)
freq_high=10/24; %[d] 10-hour cycle --> (2.4 cyc/day)
t_win=29.5307; %[d]   window size
n_shift=0.2; %0.1; % percent of each shift of data segment (0.2 = 80% overlap b/w segements)
per_gap=0.2; %0.1; % percent gap of the time window (﻿If the gap comprises more than XX% of the segment, no further analysis is attempted for the segment")

%PlotWellResponse.m
PlotTime3=PlotTime1; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime4=PlotTime2; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of the dataset timeframe.
    
%avg_amp_phase.m
PlotTime5='2012-10-01'; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime6='2012-11-15'; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of
    % the timeframe during which the user would like to view the average response.