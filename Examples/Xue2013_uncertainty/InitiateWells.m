% initialize wells%%    
% Build the Well Structure
%   - The user will edit the inputs for each well in the dataset to build the Well Structure for the test.
%     If the water level data for all wells are formatted the same (e.g.,
%     format 'S), simply copy that section and edit the iWell number and file
%     names for each well.
%   - Comment out iWell blocks for formats you aren't using.

% Template for wells with format 'S'
iWell=1;
        WellInitial(iWell).filename='data/processed/waterlevel_preprocessed_hackfix_SERIAL.csv';%'waterlevel_preprocessed_hackfix-abridged.csv';%'waterlevel_preprocessed_hackfix.csv'; % name of data file for Well 1
        WellInitial(iWell).Wellname='Well_WFSD-1'; % name of Well 1
        WellInitial(iWell).format='S'; % format code for Well 1
        WellInitial(iWell).startrow=152;%610000;%152; % (for some reason, doesn't seem to like a startrow of 1 or 2) start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=1e9;%inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='wfsd-1_tides.txt'; % name of synthetic tide for Well 1
        WellInitial(iWell).nports=1; % number of ports in Well 1
        
%         % If there is more than one noise field in the data add start and
%         % end dates as necessary. Remember to update BuildWell as well.
%         % (HAVE NOT MODIFIED NOISE YET) !!!!!
        WellInitial(iWell).noisestart1='2011-07-01'; % start of noise in data
        WellInitial(iWell).noiseend1='2011-08-07'; % end of noise in data
        
        doi=2624.67; % depth to open interval in ft (800 m)
        WellInitial(iWell).dom=(doi)*0.3048; % convert depth to open (do) interval to meters