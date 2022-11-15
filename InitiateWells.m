% initialize wells%%    
% Build the Well Structure
% The user will edit the inputs for each well in the dataset to build the Well Structure for the test.


iWell=1;
        WellInitial(iWell).filename='WellS.csv'; % name of data file for Well A
        WellInitial(iWell).Wellname='Well S'; % name of Well A
        WellInitial(iWell).format='S'; % format code for Well A
        WellInitial(iWell).startrow=2; % start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='tidenameS.dat'; % name of synthetic tide for Well A
        WellInitial(iWell).nports=1; % number of ports in Well A
        
        % If there is more than one noise field in the data add start and
        % end dates as necessary. Remember to update BuildWell as well.
        WellInitial(iWell).noisestart1='2012-10-18 15:00'; % start of noise in data
        WellInitial(iWell).noiseend1='2012-10-18 15:30'; % end of noise in data
        
        do=0; % depth to open interval in ft
        WellInitial(iWell).dom=(do)*0.3048; % convert to meters

        
   iWell=2;
        WellInitial(iWell).filename='WellA.csv'; % name of data file for Well A
        WellInitial(iWell).Wellname='Well A'; % name of Well A
        WellInitial(iWell).format='A'; % format code for Well A
        WellInitial(iWell).startrow=55; % start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='tidenameA.dat'; % name of synthetic tide for Well A
        WellInitial(iWell).nports=1; % number of ports in Well A
        
        % If there is more than one noise field in the data add start and
        % end dates as necessary. Remember to update BuildWell as well.
        WellInitial(iWell).noisestart1='2012-10-18 15:00'; % start of noise in data
        WellInitial(iWell).noiseend1='2012-10-18 15:30'; % end of noise in data
        
        do=0; % depth to open interval in ft
        WellInitial(iWell).dom=(do)*0.3048; % convert to meters
       
        
        
   iWell=3;
        WellInitial(iWell).filename='WellB.csv'; % name of data file for Well B
        WellInitial(iWell).Wellname='Well B'; % name of Well B
        WellInitial(iWell).format='B'; % format code for Well B
        WellInitial(iWell).startrow=22;% start row of data to be analyzed from csv file
        WellInitial(iWell).endrow=inf; % end row of data to be analyzed from csv file
        WellInitial(iWell).tidename='tidenameB.dat'; % name of synthetic tide for Well B
        WellInitial(iWell).nports=8; % number of ports in Well B
        
        % If there is more than one noise field per well in the data add
        % start and end dates as necessary. Remember to update BuildWell as well.
        WellInitial(iWell).noisestart_1='2013-03-24 17:00'; % start of noise in data, port a
        WellInitial(iWell).noiseend_1='2013-04-02 15:00'; % end of noise in data, port a
        WellInitial(iWell).noisestart_2='2013-03-24 17:00'; % start of noise in data, port b
        WellInitial(iWell).noiseend_2='2013-04-02 15:00'; % end of noise in data, port b
        WellInitial(iWell).noisestart_3='2013-03-24 17:00'; % start of noise in data, port c
        WellInitial(iWell).noiseend_3='2013-04-02 15:00'; % end of noise in data, port c
        WellInitial(iWell).noisestart_4='2013-03-24 17:00'; % start of noise in data, port d
        WellInitial(iWell).noiseend_4='2013-04-02 15:00'; % end of noise in data, port d
        WellInitial(iWell).noisestart_5='2013-03-24 17:00'; % start of noise in data, port e
        WellInitial(iWell).noiseend_5='2013-04-02 15:00'; % end of noise in data, port e
        WellInitial(iWell).noisestart_6='2013-03-24 17:00'; % start of noise in data, port f
        WellInitial(iWell).noiseend_6='2013-04-02 15:00'; % end of noise in data, port f
        WellInitial(iWell).noisestart_7='2013-03-24 17:00'; % start of noise in data, port g
        WellInitial(iWell).noiseend_7='2013-04-02 15:00'; % end of noise in data, port g
        WellInitial(iWell).noisestart_8='2013-03-24 17:00'; % start of noise in data, port h
        WellInitial(iWell).noiseend_8='2013-04-02 15:00'; % end of noise in data, port h
        
        do_1=185.5; % depth to open interval in ft  Ports a, b, c, d, e, f, g, h
        do_2=251.5;
        do_3=322.5;
        do_4=372.5;
        do_5=443.5;
        do_6=458.5;
        do_7=499.5;
        do_8=523.5;
        WellInitial(iWell).dom_1=(do_1)*0.3048; % convert to meters
        WellInitial(iWell).dom_2=(do_2)*0.3048;
        WellInitial(iWell).dom_3=(do_3)*0.3048;
        WellInitial(iWell).dom_4=(do_4)*0.3048;
        WellInitial(iWell).dom_5=(do_5)*0.3048;
        WellInitial(iWell).dom_6=(do_6)*0.3048;
        WellInitial(iWell).dom_7=(do_7)*0.3048;
        WellInitial(iWell).dom_8=(do_8)*0.3048;