% initialize wells%%    
% Build the Well Structure
% The user will edit the inputs for each well in the dataset to build the Well Structure for the test.


iWell=1;
        WellInitial(iWell).filename='WS-14_simple1.csv'; % name of data file for Well A
        WellInitial(iWell).Wellname='WS-14 simple'; % name of Well A
        WellInitial(iWell).format='S'; % format code for Well A
        WellInitial(iWell).startrow=152; % start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='SSFL_tides_WS-14.txt'; % name of synthetic tide for Well A
        WellInitial(iWell).nports=1; % number of ports in Well A
        
        % If there is more than one noise field in the data add start and
        % end dates as necessary. Remember to update BuildWell as well.
        WellInitial(iWell).noisestart1='2003/12/22'; % start of noise in data
        WellInitial(iWell).noiseend1='2003/12/27'; % end of noise in data                     
        
        do=40; % depth to open interval in ft
        WellInitial(iWell).dom=(do)*0.3048; % convert to meters

   iWell=2;
        WellInitial(iWell).filename='WS-14.csv'; % name of data file for Well A
        WellInitial(iWell).Wellname='WS-14'; % name of Well A
        WellInitial(iWell).format='A'; % format code for Well A
        WellInitial(iWell).startrow=190; % start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='SSFL_tides_WS-14.txt'; % name of synthetic tide for Well A
        WellInitial(iWell).nports=1; % number of ports in Well A
        
        % If there is more than one noise field in the data add start and
        % end dates as necessary. Remember to update BuildWell as well.
        WellInitial(iWell).noisestart1='2003/12/22'; % start of noise in data
        WellInitial(iWell).noiseend1='2003/12/27'; % end of noise in data                     
        
        do=40; % depth to open interval in ft
        WellInitial(iWell).dom=(do)*0.3048; % convert to meters
       
        
   iWell=3;
        WellInitial(iWell).filename='HAR-16.csv'; % name of data file for Well B
        WellInitial(iWell).Wellname='HAR-16'; % name of Well B
        WellInitial(iWell).format='B'; % format code for Well B
        WellInitial(iWell).startrow=253;% start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=33668; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='SSFL_tides_HAR-16.txt'; % name of synthetic tide for Well B
        WellInitial(iWell).nports=6; % number of ports in Well B
        
        % If there is more than one noise field per well in the data add
        % start and end dates as necessary. Remember to update BuildWell as well.
        WellInitial(iWell).noisestart_1='2003/8/25'; % start of noise in data, port 1
        WellInitial(iWell).noiseend_1='2003/5/1'; % end of noise in data, port 1
        WellInitial(iWell).noisestart_2='2003/8/25'; % start of noise in data, port 3
        WellInitial(iWell).noiseend_2='2003/5/1'; % end of noise in data, port 3
        WellInitial(iWell).noisestart_3_1='2003/10/11'; % start of noise in data, port 4                                  %%
        WellInitial(iWell).noiseend_3_1='2003/10/12 12:00:00'; % end of noise in data, port 4                             %%
        WellInitial(iWell).noisestart_3_2='2003/12/25 12:00:00'; % start of noise in data, port 4                         %%
        WellInitial(iWell).noiseend_3_2='2003/12/26 12:00:00'; % end of noise in data, port 4                             %%
        WellInitial(iWell).noisestart_3_3='2004/1/21'; % start of noise in data, port 4                                   %%
        WellInitial(iWell).noiseend_3_3='2004/1/26'; % end of noise in data, port 4                                       %%
        WellInitial(iWell).noisestart_3_4='2004/3/6'; % start of noise in data, port 4                                    %%
        WellInitial(iWell).noiseend_3_4='2004/3/11'; % end of noise in data, port 4                                       %%
        WellInitial(iWell).noisestart_3_5='2004/3/29 18:00:00'; % start of noise in data, port 4                          %%
        WellInitial(iWell).noiseend_3_5='2004/4/3'; % end of noise in data, port 4                                        %%
        WellInitial(iWell).noisestart_4_1='2003/10/21 12:00:00'; % start of noise in data, port 6                         %%
        WellInitial(iWell).noiseend_4_1='2003/10/23'; % end of noise in data, port 6                                      %%
        WellInitial(iWell).noisestart_4_2='2004/3/28'; % start of noise in data, port 6                                   %%
        WellInitial(iWell).noiseend_4_2='2004/5/1'; % end of noise in data, port 6                                        %%
        WellInitial(iWell).noisestart_5='2004/3/29 18:00:00'; % start of noise in data, port 7
        WellInitial(iWell).noiseend_5='2004/4/3'; % end of noise in data, port 7
        WellInitial(iWell).noisestart_6_1='2004/1/18'; % start of noise in data, port 8                                   %%
        WellInitial(iWell).noiseend_6_1='2004/1/23'; % end of noise in data, port 8                                       %%
        WellInitial(iWell).noisestart_6_2='2004/2/16'; % start of noise in data, port 8                                   %%
        WellInitial(iWell).noiseend_6_2='2004/2/19'; % end of noise in data, port 8                                       %%
        WellInitial(iWell).noisestart_6_3='2004/3/29 18:00:00'; % start of noise in data, port 8                          %%
        WellInitial(iWell).noiseend_6_3='2004/4/3'; % end of noise in data, port 8                                        %%
%         WellInitial(iWell).noisestart_7=''; % start of noise in data, port g
%         WellInitial(iWell).noiseend_7=''; % end of noise in data, port g
%         WellInitial(iWell).noisestart_8=''; % start of noise in data, port h
%         WellInitial(iWell).noiseend_8=''; % end of noise in data, port h
        
        do_1=2.29; % depth to open interval in ft  Ports 1, 3, 4, 6, 7, 8
        do_2=21.79;
        do_3=31.79;
        do_4=51.79;
        do_5=61.79;
        do_6=71.79;
%         do_7=;
%         do_8=;
        WellInitial(iWell).dom_1=(do_1)*0.3048; % convert to meters
        WellInitial(iWell).dom_2=(do_2)*0.3048;
        WellInitial(iWell).dom_3=(do_3)*0.3048;
        WellInitial(iWell).dom_4=(do_4)*0.3048;
        WellInitial(iWell).dom_5=(do_5)*0.3048;
        WellInitial(iWell).dom_6=(do_6)*0.3048;
%         WellInitial(iWell).dom_7=(do_7)*0.3048;
%         WellInitial(iWell).dom_8=(do_8)*0.3048;