
%Calculate Phase and Amplitude Response

##clearvars -except iWell Well
clearvars -except iWell Well t_compute  %debug
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
%[jportiz] for MC uncertainty
Ph_deg_std=[];
A_std=[];
T_std=[];


datam=(data).*0.3048; % convert to feet to meters       !!!

tUTC=Well(iWell).xi; % time for water level data, UTC
date1=Well(iWell).date1; % start time in UTC
date2=Well(iWell).date2; % end time in UTC
      
dt=Well(iWell).dt; % time interval of water data

ttide=Well(iWell).ttide; % time for synthetic tide
tide=Well(iWell).tide'; %  synthetic tide

% cut the times in the data
[dum ind_data1]=min(abs(tUTC-date1));
[dum ind_data2]=min(abs(tUTC-date2));
if (length(ind_data2)==0)
    ind_data2=length(datam);
end

%%%%the water level is not evenly sampled after noisy data is removed from raw dataset%%%%%%%%
Wcut=datam(ind_data1:ind_data2);  
TWcut=tUTC(ind_data1:ind_data2);


  tcuti=[date1:dt:date2];
  Wcuti=interp1(TWcut,Wcut,tcuti); 
##end

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
        % Determine the length of interpolated array
        len_i = length(Wcuti);
elseif length(tidecut)>length(Wcuti)
        % Interpolate the tidal data to match the size of the data being processed
        tidecut=interp1(ttide_cut,tidecut,tcuti)';
        ttide_cut=tcuti';
        % Determine the length of interpolated array
        len_i = length(tidecut)
end

%%%% MONTE-CARLO UNCERTAINTY PROPAGATION
if MCMC == 1;
##  % Pre-allocate the MC WL array 
##  mc_Wcuti = zeros(N_mc,len_i);
  % Sample WL data uncertainty from normal distribution  !!!
  % ---- Preserve values of zero if they exist
  zero_inds = find(Wcuti==0);
  wl_std = wl_precision;   % [m]  
  mc_Wcuti = Wcuti + ( randn(N_mc, length(Wcuti)) * wl_std); 
  % Re-set zero_inds back to zero
  mc_Wcuti(:,zero_inds) = 0;
##  mc_Wcuti = Wcuti + ( randn(N_mc, 1) * wl_std );
##  % Interpolate the MC generated samples
##  tcuti=[date1:dt:date2]; 
##  for i = 1:N_mc;
##    mcwci = interp1(TWcut,mc_Wcut(i,:),tcuti);  
##    mc_Wcuti(i,:) = mcwci;
##  end
end


Well(iWell).tcuti=tcuti;
Well(iWell).Wcuti=Wcuti;
% POSSIBLY ADD IN mcWcuti AT SOME POINT...      !!! !!! !!! !!! !!!
Well(iWell).tidecut=tidecut;
Well(iWell).ttide_cut=ttide_cut;
    
%%%%%%%%%%%%%%%%%%%%
% Filter and De-trend the WL and Strain Data
%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%
% Filter the MC Data
%%%%%%%%%%%%%%% 
if MCMC==1;
  display('Filtering  and detrending MC data...')    %debug
  % [[ w/in Loop ]]
  mc_Wcuti_f = zeros(N_mc, length(mc_Wcuti));  %pre-allocate
  for n = 1:N_mc;
    n
    mc_Wcuti(n,:)=(mc_Wcuti(n,:)-mean(nonzeros(mc_Wcuti(n,:)))).*(mc_Wcuti(n,:)~=0);
    mc_Wcuti(n,:)=detrend(mc_Wcuti(n,:)).*(mc_Wcuti(n,:)~=0);
    
    mc_Wcuti_f(n,:)=filtfilt(b,a,mc_Wcuti(n,:));
    mc_Wcuti_f(n,:)=mc_Wcuti_f(n,:).*(mc_Wcuti(n,:)~=0);
  end
##  % [[ Element-wise ]]
##  %   CHECK THAT THIS IS DONE CORRECTLY !!! !!! !!! !!! !!! !!! !!! !!! !!! !!!
##  mc_Wcuti=(mc_Wcuti-mean(nonzeros(mc_Wcuti))).*(mc_Wcuti~=0);
##  mc_Wcuti=detrend(mc_Wcuti).*(mc_Wcuti~=0);
##
##  mc_Wcuti_f=filtfilt(b,a,mc_Wcuti);
##  mc_Wcuti_f=mc_Wcuti_f.*(mc_Wcuti~=0);
  display('   done.')   %debug
end
%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%
% Remove Noise if specified (sets those times to 0)
%%%%%%%%%%%%%%%      
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
    %%% MC 
    if MCMC==1;
      % Zero-indices should be same for each row of MC data
      mc_Wcuti_f(:,mm1)=0;
    end
end


tidecut_f=filtfilt(b,a,tidecut);
Wcuti_f_ne0=Wcuti_f~=0;
tidecut_f=tidecut_f.*(Wcuti_f_ne0);
if MCMC ==1;
  mc_Wcuti_f_ne0=mc_Wcuti_f~=0;
end

%%%%%%%%

%%%%%%%%%%%%%%%
% [JPO] new - export the filtered, detrended data for external analysis
filename = strcat(Well(iWell).name,'_',num2str(iport),'_filtered_detrended.csv');
twocol = [ tcuti', -Wcuti_f' ];
dlmwrite(filename,twocol,'precision',20);  % csvwrite() cuts off the decimal precision
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
t_win=t_win; % days
n_shift=n_shift; % percent of each shift of data segment 
pgap=per_gap; % percent gap of the time window

%%%%%%%%%%%%%%%
% Calculate phase and amplitude responses
%%%%%%%%%%%%%%%
if MCMC==0;
  [Respon]=water_respon_time(-Wcuti_f,tidecut_f',tcuti,dt,t_win,n_shift,pgap);
  % sign convention: extension positive.
%%%%%%%%%%%%%%%
% MCMC UNCERTAINTY PROPAGATION (OPTIONAL)
%%%%%%%%%%%%%%%
elseif MCMC==1;
    display('Performing MCMC uncertainty propagation...')
    
##    % Sample WL data uncertainty from normal distribution
##    wl_std = wl_precision;   % [m]
##    mc_Wcuti_f = Wcuti_f + ( randn(N_mc, 1) * wl_std );
  for n=1:N_mc;
    fprintf('MC simulation %d of %d    (%.0f %% complete)\n', n, N_mc, n/N_mc*100); 
    [Respon]=water_respon_time(-mc_Wcuti_f(n,:),tidecut_f',tcuti,dt,t_win,n_shift,pgap);
##    [Respon]=water_respon_time(-mc_Wcuti_f,tidecut_f',tcuti,dt,t_win,n_shift,pgap);
    mcRespon(n)=Respon;
  end
    
else
  display('Valid options for MCMC are 0 or 1. Enter a valid option in UserParam.m');
end
% Save MC simulated values
##display('Saving mcRespon...')  %debug
##save mcRespon;
##load mcRespon;
##display('    done.')  %debug

if MCMC==0;
  Ph=Respon.pha_M2;
  A0=Respon.amp_M2;
  T0=Respon.time;
  Ph_deg0=Ph*180/pi;

  Ph_deg=[Ph_deg Ph_deg0];
  A=[A A0];
  T=[T T0];


%%%%%%%%%%%%%%%
% if MC simulation
%%%%%%%%%%%%%%%
else
  %Calculate the mean values for phase and amplitude response of each window
  Ph=mean([mcRespon.pha_M2],2);
  A0=mean([mcRespon.amp_M2],2);
  T0=mcRespon(1).time
##  T0=mean([mcRespon.time],2);
  Ph_deg0=Ph*180/pi;

  Ph_deg=[Ph_deg Ph_deg0];
  A=[A A0];
  T=[T T0];
  
  %Calculate the stdev for phase and amplitude response
  Ph_std=std([mcRespon.pha_M2],opt=1,dim=2);
  A0_std=std([mcRespon.amp_M2],opt=1,dim=2);
  Ph_deg0_std=Ph_std*180/pi;

  Ph_deg_std=[Ph_deg_std Ph_deg0_std];
  A_std=[A_std A0_std];
  
  %Calculate the 95% confidence intervals 
  % !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! 
  % !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! 
  % !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! 
  % !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! 
    
##  T_std=[T_std T0_std];
end

%---------------------
% TESTING
%---------------------
% plot histogram of 1st M2 pha value  m,kkl
figure;
X = [mcRespon.pha_M2](1,:) * 180/pi;
hist(X, nbins=20, norm=1.0);
##hold on;
meantext = ['mean = ' num2str(Ph_deg0(1)) '°'];
stdtext =  ['std  = ' num2str(Ph_deg0_std(1)) '°'];
text(min(X), 0.14, meantext);
text(min(X), 0.13, stdtext);
title(['MC sampled M2 phase shift [°] (Window 1) n=' num2str(N_mc)]);
##hold on;
##normX = Ph(1) + (randn(N_mc,1) *wl_std )
####normPDF = pdf('Normal', 
##normX = [min(X):0.005:max(X)];
##plot( normX, normpdf(normX, Ph(1), wl_std) )
% plot histogram of 1st M2 amp value
figure;
X = [mcRespon.amp_M2](1,:);
hist(X, nbins=20, norm=1.0);
hold on;
meantext = ['mean = ' num2str(A0(1)) ''];
stdtext =  ['std  = ' num2str(A0_std(1)) ''];
text(min(X), 0.14, meantext);
text(min(X), 0.13, stdtext);
title(['MC sampled M2 amplitude response (Window 1) n=' num2str(N_mc)])
%---------------------


% save the results in the Well structure
Well(iWell).Ph0(:,iport)=Ph_deg;
Well(iWell).A0(:,iport)=A;
Well(iWell).T0(:,iport)=T;
%(new) MC statistics
if MCMC==1;
  Well(iWell).Ph0_std(:,iport)=Ph_deg_std;
  Well(iWell).A0_std(:,iport)=A_std;
  % Save these stats to text file (mean, stdev)
  filename = strcat(Well(iWell).name,'_',num2str(iport),'_amp_stats.csv');
  twocol = [ A, A_std ];
  dlmwrite(filename,twocol,'precision',20);  % csvwrite() cuts off the decimal precision
  filename = strcat(Well(iWell).name,'_',num2str(iport),'_pha_stats.csv');
  twocol = [ Ph_deg, Ph_deg_std ];
  dlmwrite(filename,twocol,'precision',20);  % csvwrite() cuts off the decimal precision
end

end