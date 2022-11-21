% initialize wells%%    
% Build the Well Structure
% The user will edit the inputs for each well in the dataset to build the Well Structure for the test.

% Multiple wells, all formatted the same
iWell=1;
        WellInitial(iWell).filename='data/wl_DW4N.csv'; % name of data file for Well 1
        WellInitial(iWell).Wellname='Well-DW4'; % name of Well 1
        WellInitial(iWell).format='S'; % format code for Well 1
        WellInitial(iWell).startrow=1;%610000;%152; % (for some reason, doesn't seem to like a startrow of 1 or 2) start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='combinedArealStrains_dw4.txt'; % name of synthetic tide for Well 1
        WellInitial(iWell).nports=1; % number of ports in Well 1
        
%         % If there is more than one noise field in the data add start and
%         % end dates as necessary. Remember to update BuildWell as well.
%         % (HAVE NOT MODIFIED NOISE YET) !!!!!
        WellInitial(iWell).noisestart1='2011-7-1 00:00'; % start of noise in data
        WellInitial(iWell).noiseend1='2011-8-7 00:00'; % end of noise in data
        
%         do=; % depth to open interval in ft (800 m)
%         WellInitial(iWell).dom=(do)*0.3048; % convert depth to open (do) interval to meters
        WellInitial(iWell).dom=37.41+46.4; % use meters from Table1 (add WL depth + top of Open INterv. to WT) (???)
 
iWell=2;
        WellInitial(iWell).filename='data/wl_DW13.csv'; % name of data file for Well 1
        WellInitial(iWell).Wellname='Well-DW13'; % name of Well 1
        WellInitial(iWell).format='S'; % format code for Well 1
        WellInitial(iWell).startrow=1;%610000;%152; % (for some reason, doesn't seem to like a startrow of 1 or 2) start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='combinedArealStrains_dw13.txt'; % name of synthetic tide for Well 1
        WellInitial(iWell).nports=1; % number of ports in Well 1
        
%         % If there is more than one noise field in the data add start and
%         % end dates as necessary. Remember to update BuildWell as well.
%         % (HAVE NOT MODIFIED NOISE YET) !!!!!
        % From Xue2016, DW13 was jammed during August-November 2014
        WellInitial(iWell).noisestart1='2014-8-1 00:00'; % start of noise in data
        WellInitial(iWell).noiseend1='2014-12-6 00:00'; % end of noise in data
        
        do=2624.67; % depth to open interval in ft (800 m)
        WellInitial(iWell).dom=(do)*0.3048; % convert depth to open (do) interval to meters

iWell=3;
        WellInitial(iWell).filename='data/wl_DW2.csv'; % name of data file for Well 1
        WellInitial(iWell).Wellname='Well-DW2'; % name of Well 1
        WellInitial(iWell).format='S'; % format code for Well 1
        WellInitial(iWell).startrow=1;%610000;%152; % (for some reason, doesn't seem to like a startrow of 1 or 2) start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='combinedArealStrains_dw2.txt'; % name of synthetic tide for Well 1
        WellInitial(iWell).nports=1; % number of ports in Well 1
        
%         % If there is more than one noise field in the data add start and
%         % end dates as necessary. Remember to update BuildWell as well.
%         % (HAVE NOT MODIFIED NOISE YET) !!!!!
        WellInitial(iWell).noisestart1='2011-7-1 00:00'; % start of noise in data
        WellInitial(iWell).noiseend1='2011-8-7 00:00'; % end of noise in data
        
        do=2624.67; % depth to open interval in ft (800 m)
        WellInitial(iWell).dom=(do)*0.3048; % convert depth to open (do) interval to meters

iWell=4;
        WellInitial(iWell).filename='data/wl_DW11.csv'; % name of data file for Well 1
        WellInitial(iWell).Wellname='Well-DW11'; % name of Well 1
        WellInitial(iWell).format='S'; % format code for Well 1
        WellInitial(iWell).startrow=2;%610000;%152; % 
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='combinedArealStrains_dw11.txt'; % name of synthetic tide for Well 1
        WellInitial(iWell).nports=1; % number of ports in Well 1
        
%         % If there is more than one noise field in the data add start and
%         % end dates as necessary. Remember to update BuildWell as well.
%         % (HAVE NOT MODIFIED NOISE YET) !!!!!
        WellInitial(iWell).noisestart1='2011-7-1 00:00'; % start of noise in data
        WellInitial(iWell).noiseend1='2011-8-7 00:00'; % end of noise in data
        
        do=2624.67; % depth to open interval in ft (800 m)
        WellInitial(iWell).dom=(do)*0.3048; % convert depth to open (do) interval to meters
