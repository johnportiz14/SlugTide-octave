
%Calculate Phase and Amplitude Response

clearvars -except iWell Well
UserParam;

iWell
%%%%

Nports=Well(iWell).num_ports;

for iport=1:Nports,
      
    if (Nports==1)
        data=Well(iWell).yi;
    else
      data=Well(iWell).yi{iport};
    end



Ph_deg=[];
A=[];
T=[];
Well(iWell).inoise=[];


datam=(data).*0.3048; % convert to meters from feet

tUTC=Well(iWell).xi; % time for water level data, UTC
date1=Well(iWell).date1; % start time in UTC
date2=Well(iWell).date2; % end time in UTC
      
dt=Well(iWell).dt; % time interval of water data

ttide=Well(iWell).ttide; % time for synthetic tide
tide=Well(iWell).tide'; %  synthetic tide

% cut the times in the data
[dum ind_data1]=min(abs(tUTC-date1));  %[JPO] not sure where "dum" comes from originally
[dum ind_data2]=min(abs(tUTC-date2));
if (length(ind_data2)==0)
    ind_data2=length(datam);
end

%%%%the water level is not evenly sampled after noisy data is removed from raw dataset%%%%%%%%
Wcut=datam(ind_data1:ind_data2);  
TWcut=tUTC(ind_data1:ind_data2);


tcuti=[date1:dt:date2];
Wcuti=interp1(TWcut,Wcut,tcuti); 


%%%%%%%%%%%%%%%%%%%%
[dum ind_synth1]=min(abs(ttide-date1));
[dum ind_synth2]=min(abs(ttide-date2));

ttide_cut=ttide(ind_synth1:ind_synth2);
tidecut=tide(ind_synth1:ind_synth2)';


% Resize Wcuti/tcuti and tidecut to be the same size, if unequal
if length(Wcuti)>length(tidecut)
        % Interpolate the data being processed to match the size of the tidal data
        Wcuti=interp1(tcuti,Wcuti,ttide_cut);
        tcuti=ttide_cut;
elseif length(tidecut)>length(Wcuti)
        % Interpolate the tidal data to match the size of the data being processed
        tidecut=interp1(ttide_cut,tidecut,tcuti)';
        ttide_cut=tcuti';
end


Well(iWell).tcuti=tcuti;
Well(iWell).Wcuti=Wcuti;
Well(iWell).tidecut=tidecut;
Well(iWell).ttide_cut=ttide_cut;
    

%%%%%%%
Nque=1/dt/2;
f_low=1/(freq_low)/Nque; 
f_hig=1/(freq_high)/Nque; 
[b,a]=butter(2,[f_low f_hig],'bandpass');

Wcuti=(Wcuti-mean(nonzeros(Wcuti))).*(Wcuti~=0);
Wcuti=detrend(Wcuti).*(Wcuti~=0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wcuti_f=filtfilt(b,a,Wcuti);
Wcuti_f=Wcuti_f.*(Wcuti~=0);


      
    if (Nports==1)
       Well(iWell).inoise=Well(iWell).noise;
    else
       Well(iWell).inoise=Well(iWell).noise{iport};
    end

[Nnoise,dum]=size(Well(iWell).inoise);
for k=1:Nnoise,
    k1=Well(iWell).inoise(k,1);
    k2=Well(iWell).inoise(k,2);
    mm1=find(tcuti>k1 & tcuti<k2);
    Wcuti_f(mm1)=0;
end


tidecut_f=filtfilt(b,a,tidecut);
Wcuti_f_ne0=Wcuti_f~=0;
tidecut_f=tidecut_f.*(Wcuti_f_ne0);
%%%%%%%%

%%%%%%%%%%%%%%%
t_win=t_win; % days
n_shift=n_shift; % percent of each shift of data segment 
pgap=per_gap; % percent gap of the time window
[Respon]=water_respon_time(-Wcuti_f,tidecut_f',tcuti,dt,t_win,n_shift,pgap);
% sign convention: extension positive.


Ph=Respon.pha_M2;
A0=Respon.amp_M2;
T0=Respon.time;
Ph_deg0=Ph*180/pi;

Ph_deg=[Ph_deg Ph_deg0];
A=[A A0];
T=[T T0];


% save the results in the Well structure
Well(iWell).Ph0(:,iport)=Ph_deg;
Well(iWell).A0(:,iport)=A;
Well(iWell).T0(:,iport)=T;

end