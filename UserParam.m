
% Parameters used in later scripts to calculate and plot the tidal response in a set of wells.
% Inputs below are listed by the script where they are used.

%MasterWell.m
NWells=3;

%LoadTides.m
Tide_dt=1/24/60; % tidal sample spacing in units of days for the synthetic tidal model

%BuildWell.m
dt=10/24/60; % time sample spacing of the well data in units of days
tidedate='2012-10-01 00:00:00'; %date of start of synthetic tide 'YYYY-MM-DD hh:mm:ss' format
tcUTC=0; % correction from time of synthetic tide to UTC (if processed in local time).
            % may require som preprocessing if Daylight Savings Time is not already corrected for.
tide_dt=dt; % interpolate the tidal model to have the same dt as the data

%BuildWell.m, PlotOriginalWL.m, PlotWellResponse.m
cUTC=7/24; % correction from time of collected data to UTC (if collected in local time). 
            % may require som preprocessing if Daylight Savings Time is not already corrected for.
            
%PlotOriginalWL.m
PlotTime1='2012-10-01'; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime2='2013-06-30'; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of the dataset timeframe.

% %PlotOriginalWL.m, PlotWellResponse.m
% EQtime='2003-12-22 12:16:6.3'; % Date of earthquake in 'YYYY-MM-DD hh:mm:ss' format 
            % (if there is an earthquake occurence in the data).

%WaterResp_general_interp.m
freq_low=30/24; % 30 hour
freq_high=10/24; % 10 hour
t_win=29.5307; % days
n_shift=0.1; % percent of each shift of data segment 
per_gap=0.1; % percent gap of the time window

%PlotWellResponse.m
PlotTime3=PlotTime1; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime4=PlotTime2; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of the dataset timeframe.
    
%avg_amp_phase.m
PlotTime5='2012-10-01'; % Date of beginning of plot in 'YYYY-MM-DD' format
PlotTime6='2012-11-15'; % Date of ending of plot in 'YYYY-MM-DD' format
    % beginning and ending of plot should match the beginning and ending of
    % the timeframe during which the user would like to view the average response.